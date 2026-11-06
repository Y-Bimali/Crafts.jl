### Variant of the search in designabletargettree.jl that is steered by the full designability check,
### and notes, for every rule set it checks, whether the `single` and `exact` checks would have passed
### it too, to judge how well those heuristics stand in for the full check. Include only this file, and
### lib/plotting/plottargettree.jl too for `plot_heuristic_tree`.

include("../lib/designabletargettree.jl")

"""
    HeuristicNote

What `_search_designs_target_heuristics` found for one bounded rule set it checked for designability.

- `rules`, `vmf`: the rule set, and the merges that first made it from the seed.
- `parent`: the node of the search's network it was made from.
- `node`: its own node, or 0 if the full check did not pass it.
- `expanded`: whether greed kept it to have children of its own.
- `score`: minus its number of structures.
- `full`, `single`, `exact`: whether `target_isdesignable` passed it without either heuristic, with
  `single`, and with `exact`.
"""
struct HeuristicNote{R<:BindingRules}
    rules::R
    vmf::VertexMergedForm
    parent::Int
    node::Int
    expanded::Bool
    score::Int
    full::Bool
    single::Bool
    exact::Bool
end

"""
    _heuristic_checks(rules, vmf, targets, layouts; maxsize, single=missing, full=missing, counts=nothing)

Return `(full, single, exact)`, whether `target_isdesignable` passes `rules`, made from the seed by
the merges of `vmf`, without either heuristic, with `single`, and with `exact`.

The sets the three try are nested: the exact set takes one structure of each target's shape, so it is
one of the sets `single` tries, and those are the first the full check tries. A pass therefore implies
the passes of the wider checks, so the checks run from the narrowest, each only if the one before it
failed.

`single` and `full`, when not `missing`, are results already known, which are used instead of
computing them. A known full failure settles all three. `counts`, when given, is a named tuple of
`Threads.Atomic{Int}`s `exact`, `single` and `full`, each raised by one for every check computed.
"""
function _heuristic_checks(rules::BindingRules, vmf::VertexMergedForm, targets, layouts; maxsize,
                           single=missing, full=missing, counts=nothing)
    full === false && return false, false, false
    _count!(counts, :exact)
    _designable(rules, vmf, targets, layouts; maxsize, single=false, exact=true) && return true, true, true
    if ismissing(single)
        _count!(counts, :single)
        single = target_isdesignable(rules, targets; maxsize, single=true)
    end
    single && return true, true, false
    if ismissing(full)
        _count!(counts, :full)
        full = target_isdesignable(rules, targets; maxsize)
    end
    return full, false, false
end

_count!(counts, check) = isnothing(counts) || Threads.atomic_add!(counts[check], 1)

"""
    _known_full(known, m, forms, valids, scores, greed)

Return `(full, selected)` for the refinements of the `m`-th rule set expanded by the search that
`known`, the `(systems, network)` of `_search_designs_target` without heuristics, records: `selected`,
the positions in `forms` of the refinements it kept, in the order it kept them, and `full`, for each
refinement, the full check's result as far as that settles it, or `missing`.

The kept ones passed. A bounded one it did not keep failed if it would have ranked among the kept had
it passed, which is all of them when fewer than `greed` were kept, and otherwise those ranking above
the last kept, by score and then by position, as the search's stable sort ranks them. The rest are
`missing`. `forms`, `valids` and `scores` are this search's refinements of that rule set, which must
include the kept ones, bounded.
"""
function _known_full(known, m, forms, valids, scores, greed)
    systems, network = known
    children = sort!(filter(>(m), neighbors(network, m)))
    selected = map(children) do c
        f = first(_canonicalforms(systems[c]))
        i = findfirst(fs -> first(fs) == f, forms)
        isnothing(i) && error("the known search keeps a refinement of rule set $m that this one does not make")
        valids[i] || error("the known search keeps a refinement of rule set $m that is not bounded here")
        return i
    end

    full = Vector{Union{Missing,Bool}}(missing, length(forms))
    full[selected] .= true
    last = isempty(selected) ? 0 : selected[end]
    for i in eachindex(forms)
        (valids[i] && ismissing(full[i])) || continue
        if isnothing(greed) || length(selected) < greed ||
           scores[i] > scores[last] || (scores[i] == scores[last] && i < last)
            full[i] = false
        end
    end
    return full, selected
end

"""
    _search_designs_target_heuristics(seed::BindingRules, targets::AbstractVector{<:Polyform}; maxsize, maxpolys, greed, enable_rotations, verbose, isbounded_check, known=nothing, single_passes=Set(), counts=nothing)

Return `(systems, network, notes)`, searching as `_search_designs_target` with neither `single` nor
`exact`, and `notes` having a `HeuristicNote` for each bounded rule set the search checked for
designability, in the order they were checked.

Unlike `_search_designs_target`, `systems` and `network` hold every rule set the full check passed,
each child of the parent it was made from, and greed only picks which of those have children. Only
those are marked as seen, so the search expands the same rule sets `_search_designs_target` does, and
a rule set passed but not expanded under one parent can appear again under another.

The seed has to be `fully_addressable_designs` of `targets`, as for `exact`.

Checks already settled elsewhere can be skipped, giving the same result for less work:

- `known`, the `(systems, network)` of `_search_designs_target` for the same seed, targets and
  parameters, without either heuristic. The full check passed the rule sets it kept and failed those
  `_known_full` can tell it did, and the search errors if it would expand any others.
- `single_passes`, the canonical forms, `first(_canonicalforms(rules))`, of rule sets `single` passes,
  such as those kept by `_search_designs_target` with `single` or `exact` (other than the seed). The
  `single` check depends on the rule set only, not on the merges that made it, so it holds for all
  rule sets of these forms.

`counts`, when given, counts the checks computed, see `_heuristic_checks`.
"""
function _search_designs_target_heuristics(seed::BindingRules, targets::AbstractVector{<:Polyform};
                                           maxsize, maxpolys, greed, enable_rotations,
                                           verbose, isbounded_check, known=nothing, single_passes=Set(),
                                           counts=nothing)
    ns = nsites(species(seed, 1))
    all(s -> nsites(s) == ns, species(seed)) ||
        throw(ArgumentError("every species must have the same number of sites"))
    nspecies(seed) != sum(nparticles, targets) &&
        throw(ArgumentError("the seed must be the fully addressable design of the targets"))
    layouts = fully_addressable_layouts(targets)

    systems = typeof(seed)[seed]
    seen = Set(_canonicalforms(seed))
    system_network = Graph(1)
    notes = HeuristicNote{typeof(seed)}[]

    # Plain vector as a FIFO: entries are only appended and consumed in order.
    queue = Tuple{Int,typeof(seed),VertexMergedForm}[(1, seed, VertexMergedForm(nspecies(seed), species(seed, 1)))]
    head = 1

    while head <= length(queue)
        nsys = length(systems)
        n, sys, vmf = queue[head]
        # The searches expand the same rule sets in the same order, so this is the known search's node.
        m = head
        head += 1
        if verbose
            println("Evaluating system[$(nspecies(sys)), $(nbonds(sys))], nsys = $nsys, ",
                  "$(length(queue) - head + 1) in queue")
        end

        merges = _merges(nspecies(sys), ns; enable_rotations)
        new_systems = [merge_species(sys, i, j, Δψ) for (i, j, Δψ) in merges]
        new_vmfs = [merge_vertices(vmf, i, j, Δψ) for (i, j, Δψ) in merges]
        new_forms = _canonicalforms.(new_systems)
        unique_idxs = unique(i -> first(new_forms[i]), eachindex(new_forms))
        new_systems, new_vmfs, new_forms = new_systems[unique_idxs], new_vmfs[unique_idxs], new_forms[unique_idxs]

        fresh = [first(f) ∉ seen for f in new_forms]
        new_systems, new_vmfs, new_forms = new_systems[fresh], new_vmfs[fresh], new_forms[fresh]

        scores = zeros(length(new_systems))
        valids = zeros(Bool, length(new_systems))
        # `nothing` for a refinement that is not bounded, and so is never checked.
        checks = Vector{Union{NTuple{3,Bool},Nothing}}(nothing, length(new_systems))
        # :greedy hands out one refinement at a time, so the slow ones do not pile up on one thread.
        @threads :greedy for i in eachindex(new_systems)
            valids[i], scores[i] = evaluate_bounded(new_systems[i]; maxsize, maxpolys, isbounded_check)
        end
        full, selected = isnothing(known) ? (fill(missing, length(new_systems)), nothing) :
                         _known_full(known, m, new_forms, valids, scores, greed)
        @threads :greedy for i in eachindex(new_systems)
            valids[i] || continue
            single = first(new_forms[i]) in single_passes ? true : missing
            checks[i] = _heuristic_checks(new_systems[i], new_vmfs[i], targets, layouts; maxsize, single,
                                          full=full[i], counts)
            valids[i] = first(checks[i])
        end

        # Every designable refinement joins the tree, highest scoring first, and the first `greed` of
        # them are expanded. Indices are into this parent's refinements, for the notes.
        passed = findall(valids)
        ranked = passed[sortperm(scores[passed]; rev=true)]
        nexpanded = isnothing(greed) ? length(ranked) : min(greed, length(ranked))
        isnothing(selected) || ranked[1:nexpanded] == selected ||
            error("rule set $m expands other refinements than in the known search")
        node = zeros(Int, length(new_systems))
        expanded = falses(length(new_systems))
        for (j, i) in enumerate(ranked)
            node[i] = nsys + j
            push!(systems, new_systems[i])
            add_vertex!(system_network)
            add_edge!(system_network, n, nsys + j)

            j <= nexpanded || continue
            expanded[i] = true
            push!(seen, new_forms[i]...)
            push!(queue, (nsys + j, new_systems[i], new_vmfs[i]))
        end

        for (i, check) in enumerate(checks)
            isnothing(check) && continue
            push!(notes, HeuristicNote(new_systems[i], new_vmfs[i], n, node[i], expanded[i], round(Int, scores[i]),
                                       check...))
        end
    end

    verbose && println("\n\nGeneration complete: generated $(length(systems)) systems.")
    return systems, system_network, notes
end

"""
    create_designs_target_heuristics(targets::AbstractVector{<:Polyform}, particlespecies; maxsize, maxpolys, greed=nothing, enable_rotations=true, verbose=false, isbounded_check=true, known=nothing, single_passes=Set(), counts=nothing)

Return `(systems, network, notes)` from `_search_designs_target_heuristics`, seeded as
`create_designs_target` seeds its search.
"""
function create_designs_target_heuristics(targets::AbstractVector{<:Polyform}, particlespecies::ParticleSpecies;
                                          maxsize::Integer,
                                          maxpolys::Integer,
                                          greed::Union{Integer,Nothing}=nothing,
                                          enable_rotations=true,
                                          verbose=false,
                                          isbounded_check=true,
                                          known=nothing,
                                          single_passes=Set(),
                                          counts=nothing)
    seed = fully_addressable_designs(targets, particlespecies)
    return _search_designs_target_heuristics(seed, targets; maxsize, maxpolys, greed, enable_rotations, verbose,
                                             isbounded_check, known, single_passes, counts)
end

"""
    compare_trees(a::AbstractVector{<:BindingRules}, b::AbstractVector{<:BindingRules})

Return, for each number of species, how many rule sets `a` and `b`, such as the `systems` of two
searches, have, and how many of those of `a` are in `b` too, as named tuples `(nspecies, a, b, both)`
sorted by number of species. A rule set and its mirror image count as one, as the search's seen set
counts them.
"""
function compare_trees(a::AbstractVector{<:BindingRules}, b::AbstractVector{<:BindingRules})
    formsb = Set(f for r in b for f in _canonicalforms(r))
    return map(sort!(unique(nspecies.(vcat(a, b))))) do m
        la = filter(r -> nspecies(r) == m, a)
        return (; nspecies=m, a=length(la), b=count(r -> nspecies(r) == m, b),
                both=count(r -> first(_canonicalforms(r)) in formsb, la))
    end
end

"""
    sametree((systemsa, networka), (systemsb, networkb))

Return whether two searches made the same tree: the same rule sets, up to mirror images, in the same
order, with the same edges.
"""
function sametree((systemsa, networka), (systemsb, networkb))
    length(systemsa) == length(systemsb) || return false
    all(first(_canonicalforms(ra)) in _canonicalforms(rb) for (ra, rb) in zip(systemsa, systemsb)) || return false
    edgeset(g) = Set((src(e), dst(e)) for e in edges(g))
    return edgeset(networka) == edgeset(networkb)
end

"""
    regular_vmfs((systems, network); enable_rotations=true)

Return the `VertexMergedForm` of each rule set of `systems`, the result of a regular search such as
`create_designs_target`, as the search tracked it: the seed's identity form, and for every other rule
set, its parent's form merged as the first merge of the parent that makes the rule set, up to
relabeling, as the search's `unique` keeps it. `enable_rotations` must be as in the search.

Only canonical forms of merges are computed, nothing is enumerated.
"""
function regular_vmfs((systems, network); enable_rotations=true)
    seed = first(systems)
    ns = nsites(species(seed, 1))
    vmfs = Vector{VertexMergedForm}(undef, length(systems))
    vmfs[1] = VertexMergedForm(nspecies(seed), species(seed, 1))
    for c in 2:length(systems)
        # Systems are numbered as found, so the parent is the one neighbor numbered lower.
        p = minimum(neighbors(network, c))
        form = first(_canonicalforms(systems[c]))
        merges = _merges(nspecies(systems[p]), ns; enable_rotations)
        k = findfirst(m -> first(_canonicalforms(merge_species(systems[p], m...))) == form, merges)
        isnothing(k) && error("no merge of rule set $p makes rule set $c")
        vmfs[c] = merge_vertices(vmfs[p], merges[k]...)
    end
    return vmfs
end

"""
    heuristic_finds(targets, result; maxsize, exact=false, enable_rotations=true)

Return, for each rule set in `result`, the `(systems, network)` of a regular search such as
`create_designs_target` for `targets`, whether `target_isdesignable` with `single` passes it, or with
`exact`, whether the exact check passes it for the merges the search made it by, see `regular_vmfs`
(`enable_rotations` as in the search). The seed, `systems[1]`, is included.
"""
function heuristic_finds(targets::AbstractVector{<:Polyform}, result; maxsize, exact=false, enable_rotations=true)
    systems = first(result)
    vmfs = exact ? regular_vmfs(result; enable_rotations) : nothing
    layouts = exact ? fully_addressable_layouts(targets) : nothing
    finds = zeros(Bool, length(systems))
    @threads :greedy for i in eachindex(systems)
        finds[i] = exact ? _designable(systems[i], vmfs[i], targets, layouts; maxsize, single=false, exact=true) :
                   target_isdesignable(systems[i], targets; maxsize, single=true)
    end
    return finds
end

heuristic_finds(target::Polyform, result; kwargs...) = heuristic_finds([target], result; kwargs...)

"""
    heuristic_summary(notes::AbstractVector{<:HeuristicNote})

Return, for each number of species, how many rule sets of `notes` were checked, how many the full
check passed, and how many of those `single` and `exact` passed too, as named tuples sorted by number
of species.

Neither heuristic can pass a rule set the full check fails, see `_heuristic_checks`, so how many of
the full check's passes each one misses is all there is to compare.
"""
function heuristic_summary(notes::AbstractVector{<:HeuristicNote})
    return map(sort!(unique(nspecies(note.rules) for note in notes))) do m
        level = filter(note -> nspecies(note.rules) == m, notes)
        return (; nspecies=m, checked=length(level), full=count(note -> note.full, level),
                single=count(note -> note.single, level), exact=count(note -> note.exact, level))
    end
end

"""
    heuristic_tree(systems::AbstractVector{<:BindingRules}, notes::AbstractVector{<:HeuristicNote})

Return `(rules, tree, passes, occurrences)`, the tree of every rule set in `notes` the full check
passed, each once, with the seed, `systems[1]`, at vertex 1, as `_search_designs_target_heuristics`
returns `systems` and `notes`.

A rule set and its mirror image count as one, as the search's seen set counts them, and each vertex
keeps the rule set it was first found as. `tree` links each vertex to the vertex of every parent it
was found under, so a rule set found under several parents has several edges up.

`passes[v]` is the narrowest check vertex `v` passes, `:exact`, `:single` or `:subsets` (the full
check alone), as first found, or `:seed` for the seed, which the search never checks.
`occurrences[v]` lists the positions in `notes` of every time it was found. The full and `single`
checks give the same result for each, but `exact` depends on the merges that made it, so it can
differ between them.
"""
function heuristic_tree(systems::AbstractVector{<:BindingRules}, notes::AbstractVector{<:HeuristicNote})
    seed = first(systems)
    rules, passes, occurrences = [seed], [:seed], [Int[]]
    tree = Graph(1)
    vertexof = Dict(first(_canonicalforms(seed)) => 1)
    # Each node of the search's network, the seed's being 1, as its vertex of `tree`.
    vertexofnode = Dict(1 => 1)

    for (k, note) in enumerate(notes)
        note.full || continue
        g, mirror = _canonicalforms(note.rules)
        v = get(() -> get(vertexof, mirror, 0), vertexof, g)
        if v == 0
            push!(rules, note.rules)
            push!(passes, note.exact ? :exact : note.single ? :single : :subsets)
            push!(occurrences, Int[])
            add_vertex!(tree)
            v = vertexof[g] = nv(tree)
        end
        push!(occurrences[v], k)
        vertexofnode[note.node] = v
        # A parent is always expanded, and so found before its children.
        add_edge!(tree, vertexofnode[note.parent], v)
    end
    return rules, tree, passes, occurrences
end

# The values of `passes` in the order `plot_heuristic_tree` colors them, and their labels.
const HEURISTIC_PASSES = (:exact, :single, :subsets, :seed)
const HEURISTIC_PASS_LABELS = ["exact", "single", "subsets", "seed"]

# "E|S|A" for the rule sets of one level with the narrowest checks `passes`: how many pass `exact`, how
# many `single`, and how many the full check, that is all of them, the narrower passes counted in the
# wider. "seed" for the seed's level, which holds no other rule set.
function _passcounts(passes)
    checked = filter(!=(:seed), passes)
    isempty(checked) && :seed in passes && return "seed"
    nexact = count(==(:exact), checked)
    return "$nexact|$(nexact + count(==(:single), checked))|$(length(checked))"
end

# Draw the tree as `_plot_target_tree!` of plottargettree.jl does, its points colored by `passes`, but
# with each level's `_passcounts` in place of its number of rule sets, and a colorbar naming the checks.
# Returns the axis of the levels. `gl` is a `GridLayout`, left unannotated so that this file loads
# without Makie.
function _plot_heuristic_tree!(gl, data, tree::AbstractGraph, passes; xmax, title, colormap, target)
    (; ps, ns, colorvalues, colorlabel, nlevels) = data
    nclipped = count(>(xmax), ps)
    nclipped > 0 && @warn "$nclipped rule sets allow more than maxpolys = $xmax structures and lie off the plot"

    # Level axes sit at y = 1, 2, ..., the root's on top, so that y = n - nmin + 1.
    nmin, nmax = extrema(ns)
    ys = ns .- nmin .+ 1

    header = isnothing(title) ? 1 : 2
    isnothing(title) || Label(gl[1, 1:4], title; fontsize=20, font=:bold)
    # Not `L"n_{spc}"`, a macro, which would need Makie loaded before this file is.
    Label(gl[header, 1], Makie.LaTeXString(raw"$n_{spc}$"); fontsize=18)
    Label(gl[header, 3], "E|S|A")

    # The level rows are as tall as each other and without gaps, so that the labels in each one line up
    # with its level's axis in the plot spanning them. The headers sit right on top of the first one.
    levelrows = header+1:header+nlevels
    pad = 0.03xmax
    ax = Axis(gl[levelrows, 2];
              limits=(-pad, xmax + pad, 0.5, nlevels + 0.5),
              xlabel=rich("number of structures"),
              xticks=_integer_ticks(0, xmax),
              xgridcolor=(:black, 0.08), ygridvisible=false)
    hideydecorations!(ax)
    hidespines!(ax, :t, :l, :r)

    for (row, n) in zip(levelrows, nmax:-1:nmin)
        Label(gl[row, 1], string(n))
        Label(gl[row, 3], _passcounts(passes[ns .== n]))
    end
    for r in levelrows
        rowsize!(gl, r, Auto(false))
    end
    for r in first(levelrows)-1:last(levelrows)-1
        rowgap!(gl, r, 0)
    end
    isnothing(title) || rowgap!(gl, 1, TITLE_ROWGAP)
    colgap!(gl, 1, LABEL_COLGAP)
    colgap!(gl, 2, LABEL_COLGAP)

    ncat = length(HEURISTIC_PASSES)
    colors = _integer_colors(colormap, 1, ncat)
    linesegments!(ax, [Point2f(0, y) => Point2f(xmax, y) for y in 1:nlevels]; color=:gray50, linewidth=1)
    linesegments!(ax, [Point2f(ps[src(e)], ys[src(e)]) => Point2f(ps[dst(e)], ys[dst(e)]) for e in edges(tree)];
                  color=:black, linewidth=0.75)
    scatter!(ax, ps, ys; color=colorvalues, markersize=10, colors...)
    Colorbar(gl[levelrows, 4]; label=colorlabel, colormap=colors.colormap, limits=colors.colorrange,
             ticks=(1:ncat, HEURISTIC_PASS_LABELS))

    isnothing(target) || _plot_target_inset!(ax, nlevels, target)
    return ax
end

"""
    plot_heuristic_tree(rules, tree, passes; maxpolys=nothing, title="Heuristic Search Tree", filename=nothing,
                        colormap=[:seagreen, :goldenrod, :firebrick, :gray60], figsize=nothing, target=nothing)

Plot the tree `heuristic_tree` returns as `plot_target_tree` does, but with each rule set colored by
the narrowest check it passes, `passes`: `:exact`, `:single`, `:subsets` (the full check alone), or
`:seed`, in the order of `colormap`.

To the right of each level, in place of its number of rule sets, is E|S|A: how many of its rule sets
pass `exact`, how many `single`, and how many the full check, so that A is all of them. The seed's
level, which the search never checks, reads "seed".

Needs lib/plotting/plottargettree.jl included. Returns the figure, after saving it to `filename` with CairoMakie if
given.
"""
function plot_heuristic_tree(rules::AbstractVector{<:BindingRules}, tree::AbstractGraph, passes::AbstractVector{Symbol};
                             maxpolys::Union{Integer,Nothing}=nothing,
                             title="Heuristic Search Tree",
                             filename=nothing,
                             colormap=[:seagreen, :goldenrod, :firebrick, :gray60],
                             figsize=nothing,
                             target=nothing)
    length(passes) == length(rules) ||
        throw(ArgumentError("there are $(length(passes)) passes but $(length(rules)) rule sets"))
    all(in(HEURISTIC_PASSES), passes) ||
        throw(ArgumentError("passes must be among $HEURISTIC_PASSES"))

    data = _target_tree_data(rules, tree, :bonds)
    data = merge(data, (; colorvalues=[findfirst(==(p), HEURISTIC_PASSES) for p in passes], colorlabel="passes"))

    fig = Figure(; size=something(figsize, (900, _target_tree_height(data.nlevels))))
    _plot_heuristic_tree!(fig.layout, data, tree, passes; xmax=_target_tree_xmax(maxpolys, data.ps), title,
                          colormap, target)

    isnothing(filename) || save(filename, fig; backend=CairoMakie)
    return fig
end

"""
    plot_heuristic_tree((rules, tree, passes, occurrences); kwargs...)

Plot what `heuristic_tree` returns as is, see the method taking it apart.
"""
plot_heuristic_tree((rules, tree, passes, _)::Tuple; kwargs...) = plot_heuristic_tree(rules, tree, passes; kwargs...)
