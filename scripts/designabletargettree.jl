### Variant of `create_designs` in which a refinement survives the search only if some set of structures
### with the target's shape is designable, in the sense of `isdesignable`, under the rule set's own
### structure enumeration. Greed is applied after that filter, ranking the survivors by the number of
### structures they allow, as `create_designs` does.
###
### Needs Convex and Clarabel in the active environment, for example the test environment, and
### CairoMakie and GLMakie for `plot_target_tree`.

using Crafts, Roly, Graphs, Convex, Clarabel
using CairoMakie
import GLMakie
using Base.Threads

# Whichever backend loaded last is active, and the REPL shows returned figures with it. Under CairoMakie
# that opens a PNG in external viewers, so GLMakie is made active here and CairoMakie only saves files.
GLMakie.activate!()
using Crafts: generate_refinements, _canonicalforms

"""
    evaluate_bounded(rules::BindingRules; maxsize, maxpolys)

Return `(valid, score)` as `Crafts.evaluate_system` with `use_enum=true`, scoring `rules` by minus the
number of structures it allows, but abandoning the enumeration as soon as a structure of `maxsize`
particles turns up.

Reaching `maxsize` rejects the rule set anyway, so nothing past that point could change the outcome.
Prechecks for definite infinites with Roly.isbounded.
"""
function evaluate_bounded(rules::BindingRules; maxsize, maxpolys, isbounded_check=true)
    isbounded_check && (!isbounded(rules) && return (false, 0.0))
    toolarge = false
    result = polyenum(rules; maxsize, maxpolys) do _, n
        n >= maxsize || return ACCEPT
        toolarge = true
        return BREAK
    end
    valid = !toolarge && result.n < maxpolys
    return valid ? (true, convert(Float64, -result.n)) : (false, 0.0)
end

"""
    target_isdesignable(rules::BindingRules, target::Polyform; maxsize, atol=1e-6, single=false)

Return whether some nonempty set of structures of `rules`, all with the shape of `target`, is
designable.

Several structures can share the target's shape with different species assignments, and any mix of
them counts as assembling the target. The structures come from the rule set's own `AssemblySystem`,
so the indices refer to its ordering. With `single`, only sets of one structure are tried.
"""
function target_isdesignable(rules::BindingRules, target::Polyform; maxsize, atol=1e-6, single=false)
    asys = AssemblySystem(rules, EntropyModel(TreeLike()); maxsize, verbose=false)
    target_idxs = find_shape(target, structures(asys))
    n = length(target_idxs)

    for nchosen in 1:(single ? min(n, 1) : n) # size of subset, try smallest first, 1 if single, false if no target
        mask = (1 << nchosen) - 1
        while mask < 1 << n
            subset = target_idxs[[isodd(mask >> (k - 1)) for k in 1:n]]
            isdesignable(asys, subset; optimizer=Clarabel.Optimizer, atol) && return true

            # Gosper's hack: the next larger mask with the same number of set bits.
            lowbit = mask & -mask
            ripple = mask + lowbit
            mask = (((ripple ⊻ mask) >> 2) ÷ lowbit) | ripple
        end
    end
    return false
end

"""
    _search_designs_target(seed::BindingRules, target::Polyform; maxsize, maxpolys, greed, enable_rotations, verbose, single)

Return `(systems, network)` as `Crafts._search_designs`, but keeping only the refinements that are
bounded, still build `target`, and leave some set of target-shaped structures designable, and then, when `greed` is given, only the
`greed` of those allowing the fewest structures. `single` is passed on to `target_isdesignable`.
"""
function _search_designs_target(seed::BindingRules, target::Polyform;
                                maxsize, maxpolys, greed, enable_rotations,
                                verbose, single, isbounded_check)
    systems = typeof(seed)[seed]
    seen = Set(_canonicalforms(seed))
    system_network = Graph(1)

    # Plain vector as a FIFO: entries are only appended and consumed in order.
    queue = Tuple{Int,typeof(seed)}[(1, seed)]
    head = 1

    while head <= length(queue)
        nsys = length(systems)
        n, sys = queue[head]
        head += 1
        if verbose
            println("Evaluating system[$(nspecies(sys)), $(nbonds(sys))], nsys = $nsys, ",
                  "$(length(queue) - head + 1) in queue")
        end

        new_systems = generate_refinements(sys; enable_rotations)
        new_forms = _canonicalforms.(new_systems)
        unique_idxs = unique(i -> first(new_forms[i]), eachindex(new_forms))
        new_systems, new_forms = new_systems[unique_idxs], new_forms[unique_idxs]

        fresh = [first(f) ∉ seen for f in new_forms]
        new_systems, new_forms = new_systems[fresh], new_forms[fresh]

        scores = zeros(length(new_systems))
        valids = zeros(Bool, length(new_systems))
        @threads for i in eachindex(new_systems)
            valids[i], scores[i] = evaluate_bounded(new_systems[i]; maxsize, maxpolys, isbounded_check)
            # if verbose
            #     println("evaluate_bounded for $(i), $(valids[i]), $(scores[i]).")
            # end
            valids[i] && (valids[i] = target_isdesignable(new_systems[i], target; maxsize, single))
            # if verbose
            #     println("target_isdesignable for $(i), $(valids[i]).")
            # end
        end
        new_systems, new_forms, scores = new_systems[valids], new_forms[valids], scores[valids]

        ranking = sortperm(scores; rev=true)
        if !isnothing(greed)
            length(ranking) > greed && (ranking = ranking[1:greed])
        end
        new_systems, new_forms = new_systems[ranking], new_forms[ranking]

        for (j, (newsys, newform)) in enumerate(zip(new_systems, new_forms))
            push!(systems, newsys)
            push!(seen, newform...)
            push!(queue, (nsys + j, newsys))

            add_vertex!(system_network)
            add_edge!(system_network, n, nsys + j)
        end
    end

    verbose && println("\n\nGeneration complete: generated $(length(systems)) systems.")
    return systems, system_network
end

"""
    create_designs_target(target::Polyform, particlespecies; maxsize, maxpolys, greed=nothing, enable_rotations=true, verbose=false, single=false)

Return `(systems, network)` as `create_designs`, but only through rule sets under which `target` stays
designable. `greed` then keeps that many of the designable children of each rule set, preferring those
that allow the fewest structures. With `single`, designability is only checked for single target-shaped
structures, as in `target_isdesignable`.
"""
function create_designs_target(target::Polyform, particlespecies;
                               maxsize::Integer,
                               maxpolys::Integer,
                               greed::Union{Integer,Nothing}=nothing,
                               enable_rotations=true,
                               verbose=false,
                               single=false,
                               isbounded_check=true)
    seed = fully_addressable_design(target, particlespecies)
    return _search_designs_target(seed, target; maxsize, maxpolys, greed, enable_rotations, verbose, single, isbounded_check)
end

# ------ The following functions are for analysis/investigative help.

"""
    lowest_on_level(rules, n, greed)

Return the `greed` rule sets in `rules` with `n` species that allow the fewest structures,
`polyenum(rules).n`, or all of them if there are fewer, from fewest structures up.

Ties keep the order of `rules`.
"""
function lowest_on_level(rules::AbstractVector{<:BindingRules}, n::Integer, greed::Integer)
    level = filter(r -> nspecies(r) == n, rules)
    ps = zeros(Int, length(level))
    @threads for i in eachindex(level)
        ps[i] = polyenum(level[i]).n
    end
    ranking = sortperm(ps)
    return level[ranking[1:min(greed, length(ranking))]]
end

# Integer ticks from `lo` to `hi`, at most `maxticks` of them, at the smallest step of 1, 2 or 5 times a
# power of ten that allows that.
function _integer_ticks(lo::Integer, hi::Integer; maxticks=8)
    step = 1
    for scale in (10^e for e in 0:18), mantissa in (1, 2, 5)
        step = mantissa * scale
        (hi - lo) ÷ step + 1 <= maxticks && break
    end
    return collect(cld(lo, step)*step:step:hi)
end

"""
    plot_target_tree(rules, tree; colorby=:bonds, maxpolys=nothing, title="Designable Target Tree",
                     filename=nothing, colormap=:viridis, figsize=nothing)

Plot the tree of rule sets returned by `create_designs_target`, `rules` being the rule sets and `tree`
the graph on their indices.

Each level of the tree, whose rule sets share a number of species n, gets a horizontal axis, the root's
at the top, with n to its left and the number of rule sets on it to its right. Every rule set is a point
on its level's axis at the number p of structures it allows, `polyenum(rules).n`, and every edge of
`tree` a line between two points. The axes all run from 0 to `maxpolys`, by default the largest p in the
tree.

`colorby` picks what the points' colors show: `:bonds` the number of bonds, `:maxsize` the size of the
largest structure, `polyenum(rules).maxsize`, and `:multiplicity` the number of rule sets with the same n
and p, which are drawn on top of one another.

Returns the figure, after saving it to `filename` with CairoMakie if given. The REPL shows it in a
GLMakie window unless the call ends in `;`, elsewhere `display` it to get one.
"""
function plot_target_tree(rules::AbstractVector{<:BindingRules}, tree::AbstractGraph;
                          colorby::Symbol=:bonds,
                          maxpolys::Union{Integer,Nothing}=nothing,
                          title="Designable Target Tree",
                          filename=nothing,
                          colormap=:viridis,
                          figsize=nothing)
    colorby in (:bonds, :maxsize, :multiplicity) ||
        throw(ArgumentError("colorby must be :bonds, :maxsize or :multiplicity, got :$colorby"))
    nv(tree) == length(rules) ||
        throw(ArgumentError("tree has $(nv(tree)) vertices but there are $(length(rules)) rule sets"))

    ps = zeros(Int, length(rules))
    maxsizes = zeros(Int, length(rules))
    @threads for i in eachindex(rules)
        result = polyenum(rules[i])
        ps[i], maxsizes[i] = result.n, result.maxsize
    end
    ns = nspecies.(rules)

    # Level axes sit at y = 1, 2, ..., the root's on top, so that y = n - nmin + 1.
    nmin, nmax = extrema(ns)
    nlevels = nmax - nmin + 1
    ys = ns .- nmin .+ 1

    multiplicity = Dict{Tuple{Int,Int},Int}()
    foreach(np -> multiplicity[np] = get(multiplicity, np, 0) + 1, zip(ns, ps))
    colorvalues, colorlabel = if colorby == :bonds
        nbonds.(rules), L"bonds"
    elseif colorby == :maxsize
        maxsizes, "max. size"
    else
        [multiplicity[np] for np in zip(ns, ps)], "multiplicity"
    end

    xmax = something(maxpolys, max(maximum(ps), 1))
    nclipped = count(>(xmax), ps)
    nclipped > 0 && @warn "$nclipped rule sets allow more than maxpolys = $xmax structures and lie off the plot"

    figsize = something(figsize, (900, 160 + 50nlevels))
    fig = Figure(; size=figsize)
    Label(fig[1, 1:4], title; fontsize=20, font=:bold)
    Label(fig[2, 1], L"n_{spc}"; fontsize=18)
    Label(fig[2, 3], "count")

    # The level rows are as tall as each other and without gaps, so that the labels in each one line up
    # with its level's axis in the plot spanning them. The headers sit right on top of the first one.
    levelrows = 3:nlevels+2
    pad = 0.03xmax
    ax = Axis(fig[levelrows, 2];
              limits=(-pad, xmax + pad, 0.5, nlevels + 0.5),
              xlabel=rich("number of structures ", rich("p"; font=:italic)),
              xticks=_integer_ticks(0, xmax),
              xgridcolor=(:black, 0.08), ygridvisible=false)
    hideydecorations!(ax)
    hidespines!(ax, :t, :l, :r)

    for (row, n) in zip(levelrows, nmax:-1:nmin)
        Label(fig[row, 1], string(n))
        Label(fig[row, 3], string(count(==(n), ns)))
    end
    for r in levelrows
        rowsize!(fig.layout, r, Auto(false))
    end
    for r in first(levelrows)-1:last(levelrows)-1
        rowgap!(fig.layout, r, 0)
    end

    linesegments!(ax, [Point2f(0, y) => Point2f(xmax, y) for y in 1:nlevels]; color=:gray50, linewidth=1)
    linesegments!(ax, [Point2f(ps[src(e)], ys[src(e)]) => Point2f(ps[dst(e)], ys[dst(e)]) for e in edges(tree)];
                  color=:black, linewidth=0.75)

    # One color per integer value, each centered on its value in the colorbar.
    vmin, vmax = extrema(colorvalues)
    points = scatter!(ax, ps, ys; color=colorvalues, markersize=10,
                      colormap=cgrad(colormap, vmax - vmin + 1; categorical=true),
                      colorrange=(vmin - 0.5, vmax + 0.5))
    Colorbar(fig[levelrows, 4], points; label=colorlabel, ticks=_integer_ticks(vmin, vmax; maxticks=10))

    isnothing(filename) || save(filename, fig; backend=CairoMakie)
    return fig
end

"""
    plot_target_tree((rules, tree); kwargs...)

Plot the `(rules, tree)` pair returned by `create_designs_target` as is, see the method taking them apart.
"""
plot_target_tree((rules, tree)::Tuple{AbstractVector{<:BindingRules},AbstractGraph}; kwargs...) =
    plot_target_tree(rules, tree; kwargs...)
