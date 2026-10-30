### Variant of `create_designs` in which a refinement survives the search only if some set of structures
### with the target's shape is designable, in the sense of `isdesignable`, under the rule set's own
### structure enumeration. Greed is applied after that filter, ranking the survivors by the number of
### structures they allow, as `create_designs` does.
###
### Needs Convex and Clarabel in the active environment, for example the test environment.
### `plot_target_tree`, for the results, is in plottargettree.jl.

using Crafts, Roly, Graphs, Convex, Clarabel
using Base.Threads
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
    target_isdesignable(rules::BindingRules, targets::AbstractVector{<:Polyform}; maxsize, atol=1e-6, single=false)

Return whether some set of structures of `rules` with at least one structure of each target's shape,
and no others, is designable.

Several structures can share a target's shape with different species assignments, and any mix of
them counts as assembling that target. The sets tried are the unions of a nonempty subset of the
structures of each target's shape, smallest first. The structures come from the rule set's own
`AssemblySystem`, so the indices refer to its ordering. With `single`, only sets taking exactly one
structure of each target's shape are tried.

The sets with the same number of structures taken from each target's list are tried in parallel, and
once one is designable the rest are skipped and no larger sets are begun.
"""
function target_isdesignable(rules::BindingRules, targets::AbstractVector{<:Polyform}; maxsize, atol=1e-6, single=false)
    asys = AssemblySystem(rules, EntropyModel(TreeLike()); maxsize, verbose=false)
    target_idxs = [find_shape(target, structures(asys)) for target in targets]
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
    _search_designs_target(seed::BindingRules, targets::AbstractVector{<:Polyform}; maxsize, maxpolys, greed, enable_rotations, verbose, single, isbounded_check)

Return `(systems, network)` as `Crafts._search_designs`, but keeping only the refinements that are
bounded, still build every target in `targets`, and leave some set of structures with at least one of
each target's shape designable, and then, when `greed` is given, only the `greed` of those allowing
the fewest structures. `single` is passed on to `target_isdesignable`.
"""
function _search_designs_target(seed::BindingRules, targets::AbstractVector{<:Polyform};
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
        # :greedy hands out one refinement at a time, so the slow ones do not pile up on one thread.
        @threads :greedy for i in eachindex(new_systems)
            valids[i], scores[i] = evaluate_bounded(new_systems[i]; maxsize, maxpolys, isbounded_check)
            # if verbose
            #     println("evaluate_bounded for $(i), $(valids[i]), $(scores[i]).")
            # end
            valids[i] && (valids[i] = target_isdesignable(new_systems[i], targets; maxsize, single))
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
    create_designs_target(targets::AbstractVector{<:Polyform}, particlespecies; maxsize, maxpolys, greed=nothing, enable_rotations=true, verbose=false, single=false, isbounded_check=true)

Return `(systems, network)` as `create_designs`, but only through rule sets under which every target
in `targets` stays designable together, as in `target_isdesignable`. `greed` then keeps that many of
the designable children of each rule set, preferring those that allow the fewest structures. With
`single`, designability is only checked for sets of one structure per target shape.

The seed is `fully_addressable_designs(targets, particlespecies)`, so no two targets share a species.
"""
function create_designs_target(targets::AbstractVector{<:Polyform}, particlespecies::ParticleSpecies;
                               maxsize::Integer,
                               maxpolys::Integer,
                               greed::Union{Integer,Nothing}=nothing,
                               enable_rotations=true,
                               verbose=false,
                               single=false,
                               isbounded_check=true)
    seed = fully_addressable_designs(targets, particlespecies)
    return _search_designs_target(seed, targets; maxsize, maxpolys, greed, enable_rotations, verbose, single, isbounded_check)
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
