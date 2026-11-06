### Variant of the search in designabletargettree.jl that tracks each rule set as a `VertexMergedForm`
### relative to the seed, kept alongside the original for debugging and timing comparisons until it
### replaces it. Shares that file's evaluation helpers, so include only this one.
###
### The search runs a level at a time, a level being every rule set with the same number of species.
### Rule sets are grouped by canonical graph into `CanonicalRecord`s, each listing the forms found to
### describe it, so that a form met again needs no canonizing and an evaluation is never repeated.
###
### TODO, after the first draft:
### - Pruning checks each new form against every pruned form with `isdescendable`, which rebuilds the
###   new form's lookup table for every pair. Building it once per new form and reusing it across the
###   list would be cheaper, but the scan still has to stop at the first ancestor found.
### - With `reflect`, a form matched through a record's mirror key gets no map, since that map would
###   reflect every species, which `VertexAtRotation` cannot say. Such forms are still recognized
###   through `byvmf`, but are not used as sources of isomorphs. Reflection-aware maps would fix that:
###   `_canonicalforms_vmf` also returns the mirror graph's perm, a map type carries one `reflected`
###   flag for the whole map (site `k` of `y` is site `z - k` of `x`), `invert_map` keeps `z` for a
###   reflected entry, `compose_maps` reflects when exactly one of its maps does, and
###   `create_merge_isomorphs` translates a merge to `Δψ′ = zⱼ - zᵢ - Δψ` under a reflected map.
###   Worth doing only if the `:mirrored` count in `stats` turns out large.

include("../lib/designabletargettree.jl")

using Crafts: VertexMergedForm, SpeciesAtRotation, VertexAtRotation, merge_vertices, isdescendable,
              get_map_from_perm, invert_map, compose_maps, create_merge_isomorphs
using Crafts: NautyDiGraph, canonize!, labels, merge_species

"""
    BoundedStatus

Outcome of `evaluate_bounded_status`. `TooLarge` and `Unbounded` carry over to every rule set merged
from one, since a merge keeps every bond, so whatever assembled before still does. `TooManyStructures`
does not, since a merge can make distinct structures into one.
"""
@enum BoundedStatus Unevaluated Bounded TooManyStructures TooLarge Unbounded

_inherited(status::BoundedStatus) = status == TooLarge || status == Unbounded

"""
    evaluate_bounded_status(rules::BindingRules; maxsize, maxpolys, isbounded_check=true)

Return `(status, score)` as `evaluate_bounded`, but with the reason a rule set is rejected: `Unbounded`
if `Roly.isbounded` fails, `TooLarge` if a structure of `maxsize` particles turns up, and
`TooManyStructures` if the enumeration reaches `maxpolys`. The score is minus the number of
structures, or `nothing` for a rejected rule set.
"""
function evaluate_bounded_status(rules::BindingRules; maxsize, maxpolys, isbounded_check=true)
    isbounded_check && !isbounded(rules) && return Unbounded, nothing
    toolarge = false
    result = polyenum(rules; maxsize, maxpolys) do _, n
        n >= maxsize || return ACCEPT
        toolarge = true
        return BREAK
    end
    toolarge && return TooLarge, nothing
    result.n >= maxpolys && return TooManyStructures, nothing
    return Bounded, -result.n
end

"""
    _canonicalforms_vmf(rules::BindingRules)

Return `(g, mirror, perm)`: the canonical graph of `rules`, that of its mirror image, and the
permutation `canonize!` used for the first, for `get_map_from_perm`.

As `Crafts._canonicalforms`, the mirror is the transposed graph, since mirroring reverses the site
cycle of every species.
"""
function _canonicalforms_vmf(rules::BindingRules)
    g = graphrep(rules)
    mirror = NautyDiGraph(permutedims(adjacency_matrix(g)); vertex_labels=labels(g))
    perm = canonize!(g)
    canonize!(mirror)
    return g, mirror, perm
end

"""
    CanonicalRecord

Everything the search knows about one canonical graph on the current level.

- `status`, `score`: from `evaluate_bounded_status`.
- `designable`: from `target_isdesignable`, `nothing` until evaluated, and left so for a record whose
  rules are never needed.
- `vmfs`, `maps`: the forms found to describe this rule set, and for each the map taking the vertices
  of `vmfs[1]` to its own, so `maps[1]` is the identity.
- `mirrored`: forms matched only through the mirror key, with `reflect`, which have no map.
- `perm`: the `canonize!` permutation of the graph of `vmfs[1]`'s rules.
- `rules`: the rules of `vmfs[1]`, the only ones kept, since merging them is what makes the children.
- `parent`: the first parent, by position among the level's parents, to generate this.
- `partner`: without `reflect`, the record of the mirror image, which shares the result of
  `evaluate_bounded_status`.
- `evaluates_bounded`: false for a record created as a partner, so a mirror pair is evaluated once.
- `lastparent`: the last parent to list this as a child, so each parent lists it once.
- `kept`: whether greed has kept this.
"""
mutable struct CanonicalRecord{R<:BindingRules}
    status::BoundedStatus
    score::Union{Int,Nothing}
    designable::Union{Bool,Nothing}
    vmfs::Vector{VertexMergedForm}
    maps::Vector{Vector{VertexAtRotation}}
    mirrored::Vector{VertexMergedForm}
    perm::Vector{Cint}
    rules::Union{R,Nothing}
    parent::Int
    partner::Union{CanonicalRecord{R},Nothing}
    evaluates_bounded::Bool
    lastparent::Int
    kept::Bool
end
CanonicalRecord{R}(; evaluates_bounded::Bool) where {R} =
    CanonicalRecord{R}(Unevaluated, nothing, nothing, VertexMergedForm[], Vector{VertexAtRotation}[],
                       VertexMergedForm[], Cint[], nothing, 0, nothing, evaluates_bounded, 0, false)

_isvalid(r::CanonicalRecord) = r.status == Bounded && r.designable === true

_identitymap(m::Integer) = [VertexAtRotation((k, 0)) for k in 1:m]

"""
    SearchLevel

The records of one level and the two ways of finding them: `bygraph` by canonical graph, with a flag
saying whether that graph is the record's mirror form, and `byvmf` by form, which finds a form met
before without canonizing it. `records` keeps creation order, which `Dict`s do not.
"""
struct SearchLevel{R<:BindingRules,G}
    bygraph::Dict{G,Tuple{CanonicalRecord{R},Bool}}
    byvmf::Dict{VertexMergedForm,CanonicalRecord{R}}
    records::Vector{CanonicalRecord{R}}
end
SearchLevel{R,G}() where {R,G} =
    SearchLevel{R,G}(Dict{G,Tuple{CanonicalRecord{R},Bool}}(), Dict{VertexMergedForm,CanonicalRecord{R}}(),
                     CanonicalRecord{R}[])

"""
    _add_canonized!(level, vmf, rules, (g, mirror, perm), p; reflect)

Add `vmf`, just canonized from `rules` into `g`, `mirror` and `perm`, to the record of `g`, creating it
if there is none, and return `(record, map)`, the map taking the vertices of the record's first form
to `vmf`'s. The map is `nothing` if `vmf` was matched through a mirror key and so has none.

A new record gets either its mirror form as a second key, with `reflect`, or a partner record for it,
without, unless the rule set is its own mirror image.
"""
function _add_canonized!(level::SearchLevel{R}, vmf::VertexMergedForm, rules::R, (g, mirror, perm), p::Integer;
                         reflect) where {R}
    hit = get(level.bygraph, g, nothing)
    if isnothing(hit)
        r = CanonicalRecord{R}(evaluates_bounded=true)
        push!(level.records, r)
        level.bygraph[g] = (r, false)
        if mirror != g
            if reflect
                level.bygraph[mirror] = (r, true)
            else
                partner = CanonicalRecord{R}(evaluates_bounded=false)
                r.partner, partner.partner = partner, r
                push!(level.records, partner)
                level.bygraph[mirror] = (partner, false)
            end
        end
    else
        r, ismirror = hit
        if ismirror
            push!(r.mirrored, vmf)
            level.byvmf[vmf] = r
            return r, nothing
        end
    end

    level.byvmf[vmf] = r
    if isempty(r.vmfs)
        # The first form of a record, possibly a partner created empty, is the one whose rules are kept.
        r.rules, r.perm, r.parent = rules, perm, p
        vmap = _identitymap(length(vmf.vertices))
    else
        vmap = get_map_from_perm(vmf, r.vmfs[1], perm, r.perm)
    end
    push!(r.vmfs, vmf)
    push!(r.maps, vmap)
    return r, vmap
end

"""
    _add_isomorphs!(level, r, firsttochild, parent, i, j, Δψ)

Add to record `r` the merges `(i, j, Δψ)` of every isomorph of `parent`, which are isomorphic to the
child, the merge of `parent`'s first form, and return how many were new. `firsttochild` takes the
vertices of `r`'s first form to the child's.
"""
function _add_isomorphs!(level::SearchLevel, r::CanonicalRecord, firsttochild, parent::CanonicalRecord,
                         i::Integer, j::Integer, Δψ::Integer)
    length(parent.vmfs) > 1 || return 0
    n = parent.vmfs[1].nsites

    # The parent is always the first form of its record, since only that one has rules to merge, so its
    # record's maps already take the parent's vertices to its isomorphs' and need no composing. The child
    # need not be the first form of its own record, so the maps to its isomorphs are composed with the
    # one taking that record's first form to the child.
    isomorphs, isomaps = create_merge_isomorphs(parent.vmfs[1], i, j, Δψ,
                                                @view(parent.vmfs[2:end]), @view(parent.maps[2:end]))
    added = 0
    for (iso, childtoiso) in zip(isomorphs, isomaps)
        haskey(level.byvmf, iso) && continue
        level.byvmf[iso] = r
        push!(r.vmfs, iso)
        push!(r.maps, compose_maps(childtoiso, firsttochild, n))
        added += 1
    end
    return added
end

# The records in `rs` with the highest scores first, ties in their order in `rs`, keeping `greed` of
# them, or all for `nothing`.
function _best(rs::AbstractVector{<:CanonicalRecord}, greed)
    ranking = sortperm([r.score for r in rs]; rev=true)
    isnothing(greed) || (ranking = first(ranking, greed))
    return rs[ranking]
end

"""
    _select_children(children, greed, greedmode)

Return the `(record, parent)` pairs greed keeps, in the order they make up the next level, from
`children[p]`, the distinct child records of parent `p` in the order it generated them.

- `:level` keeps the `greed` highest scoring valid records of the whole level.
- `:parent` keeps the `greed` highest scoring valid children of each parent, a record kept by several
  only once.
- `:legacy` repeats `_search_designs_target`: each parent in turn keeps its `greed` highest scoring
  valid children that are not already kept and whose mirror image is not, and the next level keeps
  that order.

Under `:level` and `:parent`, a record is linked to the first parent that generated it, whichever kept
it, and the next level is sorted by score, highest first, ties in order of generation. Under `:legacy`
it is linked to the parent that kept it.
"""
function _select_children(children::AbstractVector{<:AbstractVector{R}}, greed, greedmode::Symbol) where {R<:CanonicalRecord}
    if greedmode == :legacy
        selected = Tuple{R,Int}[]
        for (p, rs) in enumerate(children)
            fresh = filter(r -> _isvalid(r) && !r.kept && !(!isnothing(r.partner) && r.partner.kept), rs)
            for r in _best(fresh, greed)
                r.kept = true
                push!(selected, (r, p))
            end
        end
        return selected
    end

    if greedmode == :level
        kept = _best(filter(_isvalid, unique(reduce(vcat, children))), greed)
    else
        kept = R[]
        for rs in children, r in _best(filter(_isvalid, rs), greed)
            r.kept || push!(kept, r)
            r.kept = true
        end
        kept = _best(kept, nothing)
    end
    foreach(r -> r.kept = true, kept)
    return [(r, r.parent) for r in kept]
end

"""
    _search_designs_target_vmf(seed::BindingRules, targets::AbstractVector{<:Polyform}; maxsize, maxpolys, greed, greedmode, enable_rotations, verbose, single, isbounded_check, reflect=false, prune=true, exact=false, stats=nothing)

Return `(systems, network)` as `_search_designs_target`, searching a level at a time.

Each level runs in phases:
- A, serially: every merge of every parent is made as a form. A form met before on this level is
  looked up in `byvmf`. Otherwise, it is dropped if it can be reached by merges from a form already
  rejected for an inherited reason or as not designable (with `prune`), and is otherwise built as rules,
  canonized, and filed under its record, together with the same merge of each isomorph of its parent.
- B, threaded: `evaluate_bounded_status` once per mirror pair, the partner sharing the result.
- D, threaded: `target_isdesignable` for every bounded record with rules.
Then the rejected records' forms are added to the pruned forms, and greed picks the next level, see
`_select_children` for `greedmode`, which is `:level`, `:parent` or `:legacy`.

With `reflect`, a rule set and its mirror image are one record, and a structure counts as having a
target's shape if its mirror image does. With `exact`, a record's designability is checked for the
structures its first form implies for the targets, see `implied_structure`, which takes the seed to be
`fully_addressable_designs` of `targets`, and `reflect` must be false. Every species must have the
same number of sites. A `Dict{Symbol,Real}` passed as `stats` is added to with counts and times, per
phase.
"""
function _search_designs_target_vmf(seed::BindingRules, targets::AbstractVector{<:Polyform};
                                    maxsize, maxpolys, greed, greedmode::Symbol, enable_rotations,
                                    verbose, single, isbounded_check, reflect=false, prune=true, exact=false, stats=nothing)
    greedmode in (:level, :parent, :legacy) ||
        throw(ArgumentError("greedmode must be :level, :parent or :legacy, got :$greedmode"))
    n = nsites(species(seed, 1))
    all(s -> nsites(s) == n, species(seed)) ||
        throw(ArgumentError("every species must have the same number of sites"))
    exact && reflect && throw(ArgumentError("reflect must be false when exact is true"))
    exact && nspecies(seed) != sum(nparticles, targets) &&
        throw(ArgumentError("with exact, the seed must be the fully addressable design of the targets"))
    layouts = exact ? fully_addressable_layouts(targets) : nothing

    R = typeof(seed)
    G = typeof(first(_canonicalforms_vmf(seed)))

    seedrecord = CanonicalRecord{R}(evaluates_bounded=false)
    push!(seedrecord.vmfs, VertexMergedForm(nspecies(seed), species(seed, 1)))
    push!(seedrecord.maps, _identitymap(nspecies(seed)))
    seedrecord.rules = seed

    systems = R[seed]
    system_network = Graph(1)
    parents = [seedrecord]
    parentnodes = [1]
    pruned = VertexMergedForm[]

    while !isempty(parents)
        nchildspecies = length(first(parents).vmfs[1].vertices) - 1
        level = SearchLevel{R,G}()
        children = [CanonicalRecord{R}[] for _ in parents]
        counts = Dict{Symbol,Real}()
        count!(key, x=1) = (counts[key] = get(counts, key, 0) + x)

        counts[:time_merging] = @elapsed for (p, parent) in enumerate(parents)
            vmf, rules = parent.vmfs[1], parent.rules
            for (i, j, Δψ) in _merges(length(vmf.vertices), n; enable_rotations)
                child = merge_vertices(vmf, i, j, Δψ)
                count!(:merges)
                r = get(level.byvmf, child, nothing)
                if isnothing(r)
                    if prune && any(bad -> isdescendable(child, bad), pruned)
                        count!(:pruned)
                        continue
                    end
                    childrules = merge_species(rules, i, j, Δψ)
                    r, firsttochild = _add_canonized!(level, child, childrules, _canonicalforms_vmf(childrules), p;
                                                      reflect)
                    count!(:canonized)
                    if isnothing(firsttochild)
                        count!(:mirrored)
                    else
                        count!(:isomorphs, _add_isomorphs!(level, r, firsttochild, parent, i, j, Δψ))
                    end
                else
                    count!(:known)
                end

                r.lastparent == p && continue
                r.lastparent = p
                push!(children[p], r)
            end
        end

        bounding = filter(r -> r.evaluates_bounded, level.records)
        counts[:evaluate_bounded] = length(bounding)
        # :greedy hands out one record at a time, so the slow ones do not pile up on one thread.
        counts[:time_bounded] = @elapsed @threads :greedy for r in bounding
            r.status, r.score = evaluate_bounded_status(r.rules; maxsize, maxpolys, isbounded_check)
            # A mirror image allows the same structures, reflected, so shares the result.
            isnothing(r.partner) || ((r.partner.status, r.partner.score) = (r.status, r.score))
        end

        designing = filter(r -> !isnothing(r.rules) && r.status == Bounded, level.records)
        counts[:target_isdesignable] = length(designing)
        counts[:time_designable] = @elapsed @threads :greedy for r in designing
            # `r.rules` are those of `vmfs[1]`, so with `exact` the implied structures come from that form.
            r.designable = _designable(r.rules, r.vmfs[1], targets, layouts; maxsize, single, exact, reflect)
        end

        # Every form of a record is pruned on the result for `vmfs[1]`.
        if prune
            for r in level.records
                _inherited(r.status) || r.designable === false || continue
                append!(pruned, r.vmfs)
                append!(pruned, r.mirrored)
            end
        end

        selected = _select_children(children, greed, greedmode)
        nsys = length(systems)
        for (k, (r, p)) in enumerate(selected)
            push!(systems, r.rules)
            add_vertex!(system_network)
            add_edge!(system_network, parentnodes[p], nsys + k)
        end
        parents = first.(selected)
        parentnodes = collect(nsys .+ eachindex(selected))

        counts[:records] = length(level.records)
        counts[:kept] = length(selected)
        isnothing(stats) || mergewith!(+, stats, counts)
        if verbose
            summary = ("$k = $(v isa AbstractFloat ? round(v; digits=2) : v)" for (k, v) in sort!(collect(counts); by=first))
            println("$nchildspecies species: ", join(summary, ", "))
        end
    end

    verbose && println("\nGeneration complete: generated $(length(systems)) systems.")
    return systems, system_network
end

"""
    create_designs_target_vmf(targets::AbstractVector{<:Polyform}, particlespecies; maxsize, maxpolys, greed=nothing, greedmode, enable_rotations=true, verbose=false, single=false, isbounded_check=true, reflect=false, prune=true, exact=false, stats=nothing)

Return `(systems, network)` as `create_designs_target`, but from `_search_designs_target_vmf`, which
takes `greedmode`, `reflect`, `prune` and `stats`.
"""
function create_designs_target_vmf(targets::AbstractVector{<:Polyform}, particlespecies::ParticleSpecies;
                                   maxsize::Integer,
                                   maxpolys::Integer,
                                   greed::Union{Integer,Nothing}=nothing,
                                   greedmode::Symbol,
                                   enable_rotations::Bool=true,
                                   verbose::Bool=false,
                                   single::Bool=false,
                                   isbounded_check::Bool=true,
                                   reflect::Bool=false,
                                   prune::Bool=true,
                                   exact::Bool=false,
                                   stats=nothing)
    seed = fully_addressable_designs(targets, particlespecies)
    return _search_designs_target_vmf(seed, targets; maxsize, maxpolys, greed, greedmode, enable_rotations, verbose,
                                      single, isbounded_check, reflect, prune, exact, stats)
end

"""
    create_designs_target_vmf(target::Polyform, particlespecies; kwargs...)

Search for the single target `target`, see the method taking a vector of targets.
"""
create_designs_target_vmf(target::Polyform, particlespecies::ParticleSpecies; kwargs...) =
    create_designs_target_vmf([target], particlespecies; kwargs...)
