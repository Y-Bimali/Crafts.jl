### Variant of `create_designs` in which a refinement survives the search only if some set of structures
### with the target's shape is designable, in the sense of `isdesignable`, under the rule set's own
### structure enumeration. Greed is applied after that filter, ranking the survivors by the number of
### structures they allow, as `create_designs` does.
###
### Needs Convex and Clarabel in the active environment, as the scripts environment has them, run with
### `julia --project=scripts`. `plot_target_tree`, for the results, is in plotting/plottargettree.jl.

using Crafts, Roly, Graphs, Convex, Clarabel
using Base.Threads
using Crafts: merge_species, _canonicalforms, canonical_id
using Crafts: VertexMergedForm, merge_vertices, placements

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

# All subsets of `v` with `k` elements, in the order of their bit masks.
function _subsets(v::AbstractVector, k::Integer)
    n = length(v)
    subsets = Vector{eltype(v)}[]
    mask = (1 << k) - 1
    while mask < 1 << n
        push!(subsets, v[[isodd(mask >> (j - 1)) for j in 1:n]])

        # Gosper's hack: the next larger mask with the same number of set bits.
        lowbit = mask & -mask
        ripple = mask + lowbit
        mask = (((ripple ⊻ mask) >> 2) ÷ lowbit) | ripple
    end
    return subsets
end

"""
    reflect_polyform(p::Polyform)

Return the mirror image of the 2D polyform `p`, meaningful only as a shape, as `find_shape` reads it.

Reflecting a particle reverses the order of its sites, taken here as site `k` of `n` becoming site
`mod1(2 - k, n)`. That is the same particle relabeled only for a species that is mirror symmetric
under this relabeling, as regular polygons such as `UnitSquare` and `UnitTriangle` are, which this
assumes. The mirror is rebuilt bond by bond, breadth first from the first particle, under rules in
which every particle is its own species and only the reflected bonds are allowed, as
`fully_addressable_design` makes them.
"""
function reflect_polyform(p::Polyform)
    dimension(p) == 2 || throw(ArgumentError("only 2D polyforms can be reflected, got $(dimension(p))D"))
    # A lone regular polygon is its own mirror image, and has no bonds to build rules from.
    nparticles(p) == 1 && return p

    spcs = [species(bindingrules(p), part.speciesindex) for part in p.particles]
    mirrored(particle, site) = mod1(2 - site, nsites(spcs[particle]))
    bondlist = map(collect(bonds(p))) do bnd
        lhs, rhs = first(bnd), last(bnd)
        (lhs.particle, mirrored(lhs.particle, lhs.site), rhs.particle, mirrored(rhs.particle, rhs.site))
    end
    rules = BindingRules(permutedims(reduce(hcat, collect.(bondlist))), spcs)

    neighbors = [Tuple{Int,Int,Int}[] for _ in 1:nparticles(p)]
    for (a, sa, b, sb) in bondlist
        push!(neighbors[a], (sa, b, sb))
        push!(neighbors[b], (sb, a, sa))
    end

    # Particles of the mirror are numbered in the order they are attached, which `index` records.
    q = Roly.Polyform(rules, 1)
    index = Dict(1 => 1)
    queue = [1]
    head = 1
    while head <= length(queue)
        a = queue[head]
        head += 1
        for (sa, b, sb) in neighbors[a]
            haskey(index, b) && continue
            attached = Roly.raise!(q, bindingsite(q, ParticleSiteLoc(index[a], sa)), SpeciesSiteLoc(b, sb))
            ismissing(attached) && error("could not attach particle $b of the mirror image")
            index[b] = nparticles(q)
            push!(queue, b)
        end
    end
    return q
end

"""
    structureid(p::Polyform)

Return the canonical id of `p` with its species labels kept, so that two polyforms share an id
exactly when they are the same structure, down to which species sit where and how each is turned.
"""
structureid(p::Polyform) = canonical_id(graphrep(p))

"""
    find_structure(p::Polyform, strs::AbstractVector{<:Polyform})

Return the positions in `strs` of the structures that are `p`, species and their rotations included,
as `find_shape` does for shapes.

`strs` is assumed sorted by number of particles, as `polygen` returns it, so the scan stops once the
structures grow past `p`.
"""
function find_structure(p::Polyform, strs::AbstractVector{<:Polyform})
    needle = structureid(p)
    n = nparticles(p)

    matches = Int[]
    for (i, str) in enumerate(strs)
        m = nparticles(str)
        m > n && break
        m < n && continue
        structureid(str) == needle && push!(matches, i)
    end
    return matches
end

"""
    target_isdesignable(rules::BindingRules, targets::AbstractVector{<:Polyform}; maxsize, atol=1e-6, single=false, reflect=false, exact=false)

Return whether some set of structures of `rules` with at least one structure of each target's shape,
and no others, is designable.

Several structures can share a target's shape with different species assignments, and any mix of
them counts as assembling that target. The sets tried are the unions of a nonempty subset of the
structures of each target's shape, smallest first. The structures come from the rule set's own
`AssemblySystem`, so the indices refer to its ordering. With `single`, only sets taking exactly one
structure of each target's shape are tried.

With `reflect`, a structure counts as having a target's shape if it has the shape of the target or of
its mirror image, as when assembled 2D structures can flip over in 3D. See `reflect_polyform` for what
that assumes of the species.

With `exact`, the targets are structures of `rules` themselves, as `implied_structure` builds them,
and only the set of exactly those is tried. Each has to be among the structures of `rules`, which
assembles it whenever it has fewer than `maxsize` particles, so a missing one throws. `reflect` must
then be false, and `single` changes nothing, since the set already takes one structure per target.

The sets with the same number of structures taken from each target's list are tried in parallel, and
once one is designable the rest are skipped and no larger sets are begun.
"""
function target_isdesignable(rules::BindingRules, targets::AbstractVector{<:Polyform}; maxsize, atol=1e-6,
                             single=false, reflect=false, exact=false)
    exact && reflect && throw(ArgumentError("reflect must be false when exact is true"))
    asys = AssemblySystem(rules, EntropyModel(TreeLike()); maxsize, verbose=false)
    if exact
        target_idxs = [find_structure(target, structures(asys)) for target in targets]
        all(idxs -> length(idxs) == 1, target_idxs) ||
            error("each exact target should be exactly one structure of the rules, but they are found $(length.(target_idxs)) times")
        return isdesignable(compositionmatrix(asys), unique(reduce(vcat, target_idxs)); optimizer=Clarabel.Optimizer, atol)
    end

    target_idxs = [find_shape(target, structures(asys)) for target in targets]
    if reflect
        reflected_target_idxs = [find_shape(reflect_polyform(target), structures(asys)) for target in targets]
        target_idxs = [sort!(union(idxs, ridxs)) for (idxs, ridxs) in zip(target_idxs, reflected_target_idxs)]
    end
    any(isempty, target_idxs) && return false

    # Computed before the threads start, so that they share only this matrix, which none of them writes.
    M = compositionmatrix(asys)
    # Qualified, since Makie, loaded by plottargettree.jl, exports an `Atomic` of its own.
    found = Threads.Atomic{Bool}(false)

    # Number of structures taken from each target's list, with the smallest totals first.
    counts = Iterators.product((1:(single ? 1 : length(idxs)) for idxs in target_idxs)...)
    for nchosen in sort!(vec(collect(counts)); by=sum)
        parts = [_subsets(idxs, k) for (idxs, k) in zip(target_idxs, nchosen)]
        @threads :greedy for choice in Iterators.product(parts...)
            found[] && continue
            subset = unique(reduce(vcat, choice))
            isdesignable(M, subset; optimizer=Clarabel.Optimizer, atol) && (found[] = true)
        end
        found[] && return true
    end
    return false
end

"""
    _merges(m, n; enable_rotations)

Return every merge `(i, j, Δψ)` of `m` species with `n` sites each, in the order `generate_refinements`
makes them.
"""
_merges(m::Integer, n::Integer; enable_rotations) =
    [(i, j, Δψ) for i in 1:m for j in i+1:m for Δψ in (enable_rotations ? (0:n-1) : (0:0))]

"""
    _designable(rules::BindingRules, vmf::VertexMergedForm, targets, layouts; maxsize, single, exact, reflect=false)

Return `target_isdesignable` for `rules`, made from the seed by the merges of `vmf`: of the
structures `vmf` implies for `layouts` with `exact`, see `implied_structure`, and of the shapes of
`targets` without, with `single` passed on.
"""
function _designable(rules::BindingRules, vmf::VertexMergedForm, targets, layouts; maxsize, single, exact,
                     reflect=false)
    exact || return target_isdesignable(rules, targets; maxsize, single, reflect)
    implied = [implied_structure(rules, vmf, layout) for layout in layouts]
    return target_isdesignable(rules, implied; maxsize, exact)
end

"""
    _search_designs_target(seed::BindingRules, targets::AbstractVector{<:Polyform}; maxsize, maxpolys, greed, enable_rotations, verbose, single, isbounded_check, exact=false)

Return `(systems, network)` as `Crafts._search_designs`, but keeping only the refinements that are
bounded, still build every target in `targets`, and leave some set of structures with at least one of
each target's shape designable, and then, when `greed` is given, only the `greed` of those allowing
the fewest structures. `single` is passed on to `target_isdesignable`.

Each rule set is tracked with the `VertexMergedForm` of the merges that first made it from the seed.
With `exact`, the set checked for designability is instead the structures that form implies for the
targets, see `implied_structure`, which takes the seed to be `fully_addressable_designs` of `targets`.
Every species must have the same number of sites.
"""
function _search_designs_target(seed::BindingRules, targets::AbstractVector{<:Polyform};
                                maxsize, maxpolys, greed, enable_rotations,
                                verbose, single, isbounded_check, exact=false)
    ns = nsites(species(seed, 1))
    all(s -> nsites(s) == ns, species(seed)) ||
        throw(ArgumentError("every species must have the same number of sites"))
    exact && nspecies(seed) != sum(nparticles, targets) &&
        throw(ArgumentError("with exact, the seed must be the fully addressable design of the targets"))
    layouts = exact ? fully_addressable_layouts(targets) : nothing

    systems = typeof(seed)[seed]
    seen = Set(_canonicalforms(seed))
    system_network = Graph(1)

    # Plain vector as a FIFO: entries are only appended and consumed in order.
    queue = Tuple{Int,typeof(seed),VertexMergedForm}[(1, seed, VertexMergedForm(nspecies(seed), species(seed, 1)))]
    head = 1

    while head <= length(queue)
        nsys = length(systems)
        n, sys, vmf = queue[head]
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
        # :greedy hands out one refinement at a time, so the slow ones do not pile up on one thread.
        @threads :greedy for i in eachindex(new_systems)
            valids[i], scores[i] = evaluate_bounded(new_systems[i]; maxsize, maxpolys, isbounded_check)
            # if verbose
            #     println("evaluate_bounded for $(i), $(valids[i]), $(scores[i]).")
            # end
            valids[i] && (valids[i] = _designable(new_systems[i], new_vmfs[i], targets, layouts; maxsize, single, exact))
            # if verbose
            #     println("target_isdesignable for $(i), $(valids[i]).")
            # end
        end
        new_systems, new_vmfs, new_forms, scores = new_systems[valids], new_vmfs[valids], new_forms[valids], scores[valids]

        ranking = sortperm(scores; rev=true)
        if !isnothing(greed)
            length(ranking) > greed && (ranking = ranking[1:greed])
        end
        new_systems, new_vmfs, new_forms = new_systems[ranking], new_vmfs[ranking], new_forms[ranking]

        for (j, (newsys, newvmf, newform)) in enumerate(zip(new_systems, new_vmfs, new_forms))
            push!(systems, newsys)
            push!(seen, newform...)
            push!(queue, (nsys + j, newsys, newvmf))

            add_vertex!(system_network)
            add_edge!(system_network, n, nsys + j)
        end
    end

    verbose && println("\n\nGeneration complete: generated $(length(systems)) systems.")
    return systems, system_network
end

"""
    _search_designs_target(seed::BindingRules, target::Polyform; kwargs...)

Search for the single target `target`, see the method taking a vector of targets.
"""
_search_designs_target(seed::BindingRules, target::Polyform; kwargs...) =
    _search_designs_target(seed, [target]; kwargs...)

"""
    fully_addressable_designs(targets::AbstractVector{<:Polyform}, particlespecies::ParticleSpecies)

Return the fully addressable designs of every target in `targets`, as `fully_addressable_design`,
side by side in one rule set, with species kept apart.

The species of each target are numbered after those of the ones before it, so with 4 particles in the
first target and 5 in the second, the second's become species 5 through 9, and its bonds are
renumbered to match. No bonds are added between the targets.
"""
function fully_addressable_designs(targets::AbstractVector{<:Polyform}, particlespecies::ParticleSpecies)
    rulesets = [fully_addressable_design(target, particlespecies) for target in targets]
    bondlist = Matrix{Int}(undef, 0, 4)
    offset = 0
    for rules in rulesets
        rows = [[a.species + offset, a.site, b.species + offset, b.site]
                for (sites1, sites2) in Roly.bonded_sites(rules) for a in sites1 for b in sites2]
        isempty(rows) || (bondlist = vcat(bondlist, permutedims(reduce(hcat, rows))))
        offset += nspecies(rules)
    end
    return BindingRules(bondlist, reduce(vcat, collect.(species.(rulesets))))
end

"""
    fully_addressable_layouts(targets::AbstractVector{<:Polyform})

Return, for each target in `targets`, how its particles became species of
`fully_addressable_designs(targets, …)`: a named tuple of `species`, particle `k` having become seed
species `species[k]`, and `bonds`, the target's bonds as rows `[particle, site, particle, site, twist]`,
as `Roly.todict` writes them.

A seed species keeps the site numbering of the particle it was made from, so the bonds also say how
the seed species fit together into the target.
"""
function fully_addressable_layouts(targets::AbstractVector{<:Polyform})
    offset = 0
    return map(targets) do target
        layout = (species=offset .+ (1:nparticles(target)), bonds=Roly.todict(target; rules=false)["bonds"])
        offset += nparticles(target)
        return layout
    end
end

"""
    implied_structure(rules::BindingRules, vmf::VertexMergedForm, layout)

Return the structure of `rules` that the target of `layout`, one of `fully_addressable_layouts`,
becomes once the seed species are merged as `vmf` says, `rules` being the rule set those merges make.

Each particle is the merged species its seed species went into, turned so that its sites bond as the
target's did, see `placements`. The particles are in the target's order. A merge keeps every bond,
so `rules` assembles this structure.
"""
function implied_structure(rules::BindingRules, vmf::VertexMergedForm, layout)
    place, n = placements(vmf), vmf.nsites
    seed = place[layout.species]
    rows = [[a, mod1(sa - seed[a].rotation, n), b, mod1(sb - seed[b].rotation, n), t]
            for (a, sa, b, sb, t) in layout.bonds]
    d = Dict{String,Any}("species" => [s.vertex for s in seed], "bonds" => rows)
    return Roly.fromdict(Polyform, d; rules)
end

"""
    create_designs_target(targets::AbstractVector{<:Polyform}, particlespecies; maxsize, maxpolys, greed=nothing, enable_rotations=true, verbose=false, single=false, isbounded_check=true, exact=false)

Return `(systems, network)` as `create_designs`, but only through rule sets under which every target
in `targets` stays designable together, as in `target_isdesignable`. `greed` then keeps that many of
the designable children of each rule set, preferring those that allow the fewest structures. With
`single`, designability is only checked for sets of one structure per target shape. With `exact`,
designability is checked for the structures the merges imply, see `_search_designs_target`.

The seed is `fully_addressable_designs(targets, particlespecies)`, so no two targets share a species.
"""
function create_designs_target(targets::AbstractVector{<:Polyform}, particlespecies::ParticleSpecies;
                               maxsize::Integer,
                               maxpolys::Integer,
                               greed::Union{Integer,Nothing}=nothing,
                               enable_rotations=true,
                               verbose=false,
                               single=false,
                               isbounded_check=true,
                               exact=false)
    seed = fully_addressable_designs(targets, particlespecies)
    return _search_designs_target(seed, targets; maxsize, maxpolys, greed, enable_rotations, verbose, single, isbounded_check, exact)
end

"""
    create_designs_target(target::Polyform, particlespecies; kwargs...)

Search for the single target `target`, see the method taking a vector of targets.
"""
create_designs_target(target::Polyform, particlespecies::ParticleSpecies; kwargs...) =
    create_designs_target([target], particlespecies; kwargs...)

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
