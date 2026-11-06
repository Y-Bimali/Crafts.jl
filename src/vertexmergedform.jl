### Representation of a `BindingRules` relative to an initial one, as the result of merging species of
### the initial rules, as `merge_species` and `generate_refinements` in bindingruledesign.jl do.

"""
    SpeciesAtRotation

Species `species` of the initial rules, at rotation `rotation` relative to the reference species of the
vertex it was merged into.
"""
const SpeciesAtRotation = @NamedTuple{species::Int, rotation::Int}

"""
    VertexAtRotation

Vertex `vertex` of one form, at rotation `rotation`, as `get_map_from_perm` returns it for a vertex of
another.
"""
const VertexAtRotation = @NamedTuple{vertex::Int, rotation::Int}

"""
    VertexMergedForm

Rule set described by which species of an initial rule set have been merged together, and at which
rotations, for species that all carry `nsites` binding sites.

Each vertex is a vector of `SpeciesAtRotation`, rotations taken relative to the vertex's reference
species and lying in `0:nsites-1`. Every species `1:N` of the initial rules appears exactly once, `N`
being the total number of entries. The entries of a vertex are sorted by species, and the first one,
the reference, always has rotation zero. The vertices are sorted by their reference species, so vertex
`k` corresponds to species `k` of the merged rules.

Forms share vertex vectors with the forms they were merged from, so a form's vertex vectors are never
modified in place: build a new form instead.
"""
struct VertexMergedForm
    vertices::Vector{Vector{SpeciesAtRotation}}
    nsites::Int

    function VertexMergedForm(vertices::Vector{Vector{SpeciesAtRotation}}, nsites::Integer)
        nsites >= 1 || throw(ArgumentError("species need at least one binding site, got $nsites"))
        all(!isempty, vertices) || throw(ArgumentError("vertices must not be empty"))
        N = sum(length, vertices; init=0)
        seen = falses(N)
        for v in vertices, e in v
            1 <= e.species <= N && !seen[e.species] ||
                throw(ArgumentError("the species must be 1:$N, each appearing exactly once"))
            seen[e.species] = true
        end
        all(v -> first(v).rotation == 0, vertices) ||
            throw(ArgumentError("the first species of every vertex must have rotation zero"))
        all(v -> all(e -> 0 <= e.rotation < nsites, v), vertices) ||
            throw(ArgumentError("rotations must lie in 0:$(nsites - 1)"))
        all(v -> issorted(v; by=e -> e.species, lt=<=), vertices) ||
            throw(ArgumentError("species within a vertex must be strictly increasing"))
        issorted(vertices; by=v -> first(v).species, lt=<=) ||
            throw(ArgumentError("vertices must be strictly increasing in their first species"))
        return new(vertices, nsites)
    end
end

# Forms are normalized on construction, so two describe the same merge exactly when their contents
# agree. Without these, `==` and `hash` would compare the vertex vectors by identity.
Base.:(==)(a::VertexMergedForm, b::VertexMergedForm) = a.nsites == b.nsites && a.vertices == b.vertices
Base.hash(vmf::VertexMergedForm, h::UInt) = hash(vmf.vertices, hash(vmf.nsites, hash(:VertexMergedForm, h)))

"""
    VertexMergedForm(vertices, nsites::Integer)

Build a form from vertices given as plain `(species, rotation)` tuples, as in
`VertexMergedForm([[(1, 0), (3, 3)], [(2, 0)]], 4)`.
"""
VertexMergedForm(vertices::AbstractVector{<:AbstractVector{<:Tuple{Integer,Integer}}}, nsites::Integer) =
    VertexMergedForm([SpeciesAtRotation.(v) for v in vertices], nsites)

"""
    VertexMergedForm(n::Integer, species::ParticleSpecies)

Return the initial form on `n` copies of `species`, in which no species have been merged.
"""
VertexMergedForm(n::Integer, species::ParticleSpecies) =
    VertexMergedForm([[SpeciesAtRotation((a, 0))] for a in 1:n], nsites(species))

"""
    merge_vertices(vmf::VertexMergedForm, i, j, Δψ)

Return the form obtained by merging vertex `j` of `vmf` into vertex `i`, rotated by `Δψ` sites.
Requires `i < j`.

Follows the convention of `merge_species`, so a species at rotation `b` within vertex `j` ends up at
rotation `b + Δψ` within the merged vertex, modulo the number of sites.
"""
function merge_vertices(vmf::VertexMergedForm, i::Integer, j::Integer, Δψ::Integer)
    i < j || throw(ArgumentError("expected i < j, got i = $i and j = $j"))
    n = vmf.nsites

    vertices = copy(vmf.vertices)
    rotated_j = [SpeciesAtRotation((e.species, mod(e.rotation + Δψ, n))) for e in vertices[j]]
    vertices[i] = sort!([vertices[i]; rotated_j]; by=e -> e.species)
    deleteat!(vertices, j)
    return VertexMergedForm(vertices, n)
end

"""
    placements(vmf::VertexMergedForm)

Return a vector of `VertexAtRotation` whose `s`th entry `(vertex=k, rotation=r)` says species `s` of
the initial rules lies in vertex `k` of `vmf` at rotation `r`.

By the convention of `merge_vertices`, site `k` of the merged species is site `mod1(k + r, n)` of
initial species `s`, so site `q` of `s` is site `mod1(q - r, n)` of the merged species, where `n` is
the number of sites per species.
"""
function placements(vmf::VertexMergedForm)
    place = Vector{VertexAtRotation}(undef, sum(length, vmf.vertices; init=0))
    for (k, v) in enumerate(vmf.vertices), e in v
        place[e.species] = VertexAtRotation((k, e.rotation))
    end
    return place
end

"""
    isdescendable(a::VertexMergedForm, b::VertexMergedForm)

Return whether `a` is obtainable from `b` by one or more merges.

This is reachability by merges alone, not ancestry in a search tree: a search that prunes or
deduplicates may never reach `a` from `b` even when this returns `true`.

That is the case exactly when `a` has fewer vertices than `b` over the same initial species, every
vertex of `b` lies within a single vertex of `a`, and the rotations `a` gives that vertex's species,
taken relative to its first species, are the ones `b` gives them. Returns as soon as one of these
fails.
"""
function isdescendable(a::VertexMergedForm, b::VertexMergedForm)
    a.nsites == b.nsites || return false
    length(a.vertices) < length(b.vertices) || return false
    N = sum(length, a.vertices; init=0)
    N == sum(length, b.vertices; init=0) || return false
    n = a.nsites

    # Where each initial species sits in `a`, so each vertex of `b` is checked in a single pass. Both
    # forms cover exactly the species 1:N, so every entry gets filled.
    avertex = Vector{Int}(undef, N)
    arotation = Vector{Int}(undef, N)
    for (k, v) in enumerate(a.vertices), e in v
        avertex[e.species] = k
        arotation[e.species] = e.rotation
    end

    for v in b.vertices
        ref = first(v).species
        k = avertex[ref]
        for e in v
            avertex[e.species] == k || return false
            mod(arotation[e.species] - arotation[ref], n) == e.rotation || return false
        end
    end
    return true
end

"""
    get_map_from_perm(vmf1::VertexMergedForm, vmf2::VertexMergedForm, perm1, perm2)

Return a vector of `VertexAtRotation` whose `x`th entry `(vertex=y, rotation=z)` means vertex `x` of
`vmf2` is vertex `y` of `vmf1` at rotation `z`, so that site `k` of `y` is site `mod1(k + z, n)` of
`x`, where `n` is the number of sites per species.

The map takes the vertices of `vmf2`, the second argument, to those of `vmf1`, the first. So the map
taking the vertices of a form `a` to those of a form `b` is `get_map_from_perm(b, a, perm_b, perm_a)`,
which is the direction `create_merge_isomorphs` expects.

`perm1` and `perm2` are the permutations `canonize!` returns for `graphrep` of the rule sets the two
forms describe, which must be isomorphic. Vertex `k` of either canonical graph is vertex `perm[k]` of
the original, so original vertices with the same canonical position correspond.

Assumes the layout `graphrep` gives a `BindingRules` with one color per site: vertex
`(s - 1) * n + k` is site `k` of species `s`, and species `s` is vertex `s` of the form.
"""
function get_map_from_perm(vmf1::VertexMergedForm, vmf2::VertexMergedForm,
                           perm1::AbstractVector{<:Integer}, perm2::AbstractVector{<:Integer})
    m = length(vmf2.vertices)
    length(vmf1.vertices) == m && vmf1.nsites == vmf2.nsites ||
        throw(ArgumentError("isomorphic forms must have the same number of vertices and sites"))
    n = vmf2.nsites

    canonpos2 = invperm(perm2)
    return map(1:m) do x
        # Site 1 of vertex x is enough, the rest of its cycle follows by rotation.
        v1 = perm1[canonpos2[(x - 1) * n + 1]]
        y, site = divrem(v1 - 1, n) .+ 1
        return VertexAtRotation((y, mod(1 - site, n)))
    end
end

"""
    invert_map(vmap::AbstractVector{VertexAtRotation}, nsites::Integer)

Return the inverse of `vmap`, a map as `get_map_from_perm` returns it, so that inverting
`get_map_from_perm(vmf1, vmf2, perm1, perm2)` gives `get_map_from_perm(vmf2, vmf1, perm2, perm1)`.

Vertex `x` being vertex `y` at rotation `z` means vertex `y` is vertex `x` at rotation `-z`, modulo
`nsites`.
"""
function invert_map(vmap::AbstractVector{VertexAtRotation}, nsites::Integer)
    inverse = Vector{VertexAtRotation}(undef, length(vmap))
    for (x, e) in enumerate(vmap)
        inverse[e.vertex] = VertexAtRotation((x, mod(-e.rotation, nsites)))
    end
    return inverse
end

"""
    compose_maps(outer::AbstractVector{VertexAtRotation}, inner::AbstractVector{VertexAtRotation}, nsites::Integer)

Return the map applying `inner` and then `outer`, both maps as `get_map_from_perm` returns them, so
that composing `get_map_from_perm(vmf1, vmf2, …)` with `get_map_from_perm(vmf2, vmf3, …)`, in that
order, gives `get_map_from_perm(vmf1, vmf3, …)`.

Rotations add: if vertex `w` is vertex `x` at rotation `u` and `x` is vertex `y` at rotation `z`, then
`w` is `y` at rotation `z + u`, modulo `nsites`.
"""
function compose_maps(outer::AbstractVector{VertexAtRotation}, inner::AbstractVector{VertexAtRotation},
                      nsites::Integer)
    return map(inner) do e
        o = outer[e.vertex]
        return VertexAtRotation((o.vertex, mod(o.rotation + e.rotation, nsites)))
    end
end

"""
    create_merge_isomorphs(vmf::VertexMergedForm, i, j, Δψ, vmfs, vmaps)

Return `(isomorphs, isomaps)`, the merge `merge_vertices(vmf, i, j, Δψ)` carried over to each form of
`vmfs`, which must be isomorphic to `vmf`.

`vmaps[k]` takes the vertices of `vmf` to those of `vmfs[k]`, as `get_map_from_perm(vmfs[k], vmf, …)`
returns it. Then `isomorphs[k]` is the merge of `vmfs[k]` isomorphic to `merge_vertices(vmf, i, j, Δψ)`,
and `isomaps[k]` takes the vertices of the latter to those of `isomorphs[k]` in the same way, so the
results can be merged again.

If `vmaps[k]` takes `i` to `yᵢ` at rotation `zᵢ` and `j` to `yⱼ` at rotation `zⱼ`, joining site `s` of `i`
to site `s + Δψ` of `j` joins site `s - zᵢ` of `yᵢ` to site `s + Δψ - zⱼ` of `yⱼ`, so `yⱼ` is merged into
`yᵢ` at rotation `Δψ + zᵢ - zⱼ`, or the other way around if `yⱼ < yᵢ`.
"""
function create_merge_isomorphs(vmf::VertexMergedForm, i::Integer, j::Integer, Δψ::Integer,
                                vmfs::AbstractVector{VertexMergedForm},
                                vmaps::AbstractVector{<:AbstractVector{VertexAtRotation}})
    i < j || throw(ArgumentError("expected i < j, got i = $i and j = $j"))
    length(vmfs) == length(vmaps) || throw(ArgumentError("expected one map per form"))
    n = vmf.nsites
    m = length(vmf.vertices)

    isomorphs = Vector{VertexMergedForm}(undef, length(vmfs))
    isomaps = Vector{Vector{VertexAtRotation}}(undef, length(vmfs))
    for (k, (iso, vmap)) in enumerate(zip(vmfs, vmaps))
        (yᵢ, zᵢ), (yⱼ, zⱼ) = vmap[i], vmap[j]
        # `merge_vertices` takes the lower vertex first, and the merged vertex keeps that one's frame.
        if yᵢ < yⱼ
            i′, j′, Δψ′, zmerged = yᵢ, yⱼ, mod(Δψ + zᵢ - zⱼ, n), zᵢ
        else
            i′, j′, Δψ′, zmerged = yⱼ, yᵢ, mod(zⱼ - zᵢ - Δψ, n), mod(zⱼ - Δψ, n)
        end
        isomorphs[k] = merge_vertices(iso, i′, j′, Δψ′)

        # Both merges drop their second vertex, shifting the ones after it down by one.
        isomaps[k] = map(1:m-1) do x
            v = x < j ? x : x + 1
            v == i && return VertexAtRotation((i′, zmerged))
            y, z = vmap[v]
            return VertexAtRotation((y < j′ ? y : y - 1, z))
        end
    end
    return isomorphs, isomaps
end
