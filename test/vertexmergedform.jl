@testset "vertexmergedform" begin
    VMF = Crafts.VertexMergedForm
    named(vertices) = [Crafts.SpeciesAtRotation.(v) for v in vertices]
    plain(vmf) = [Tuple.(v) for v in vmf.vertices]

    @testset "construction" begin
        vmf = VMF([[(1, 0), (5, 3)], [(2, 0), (3, 2)], [(4, 0)]], 4)
        @test vmf.nsites == 4
        @test vmf.vertices[1][2] == (species=5, rotation=3)

        @test plain(VMF(3, UnitSquare)) == [[(1, 0)], [(2, 0)], [(3, 0)]]
        @test VMF(3, UnitTriangle).nsites == 3

        @test_throws ArgumentError VMF([[(1, 0)]], 0)
        @test_throws ArgumentError VMF([[(1, 0)], Tuple{Int,Int}[]], 4)
        @test_throws ArgumentError VMF([[(1, 1)]], 4)
        @test_throws ArgumentError VMF([[(1, 0), (2, 4)]], 4)
        @test_throws ArgumentError VMF([[(1, 0), (2, -1)]], 4)
        @test_throws ArgumentError VMF([[(2, 0), (1, 1)]], 4)
        @test_throws ArgumentError VMF([[(1, 0), (1, 1)]], 4)
        @test_throws ArgumentError VMF([[(2, 0)], [(1, 0)]], 4)
        # The species have to be exactly 1:N.
        @test_throws ArgumentError VMF([[(1, 0)], [(3, 0)]], 4)
        @test_throws ArgumentError VMF([[(0, 0)], [(1, 0)]], 4)
        @test_throws ArgumentError VMF([[(1, 0), (2, 1)], [(4, 0)]], 4)
        @test isempty(VMF(Vector{Tuple{Int,Int}}[], 4).vertices)
    end

    @testset "equality and hashing" begin
        a = VMF([[(1, 0), (3, 3)], [(2, 0)]], 4)
        b = VMF([[(1, 0), (3, 3)], [(2, 0)]], 4)
        @test a == b
        @test hash(a) == hash(b)
        @test a != VMF([[(1, 0), (3, 2)], [(2, 0)]], 4)
        @test a != VMF([[(1, 0), (3, 3)], [(2, 0)]], 5)
        @test VMF(3, UnitSquare) != VMF(3, UnitTriangle)

        # The same form reached by merging in a different order.
        initial = VMF(3, UnitSquare)
        viafirst = Crafts.merge_vertices(Crafts.merge_vertices(initial, 1, 2, 1), 1, 2, 2)
        vialast = Crafts.merge_vertices(Crafts.merge_vertices(initial, 2, 3, 1), 1, 2, 1)
        @test viafirst == vialast
        @test hash(viafirst) == hash(vialast)

        # Usable as keys, found by contents rather than by identity.
        @test Dict(a => 1)[b] == 1
        @test length(Set([a, b, viafirst, vialast])) == 2
    end

    @testset "merge_vertices" begin
        vmf = VMF([[(1, 0), (5, 3)], [(2, 0), (3, 2)], [(4, 0)]], 4)
        merged = Crafts.merge_vertices(vmf, 1, 2, 2)
        @test merged.vertices == named([[(1, 0), (2, 2), (3, 0), (5, 3)], [(4, 0)]])
        @test merged.nsites == 4

        # Vertices after the dropped one move up, keeping their order.
        @test plain(Crafts.merge_vertices(VMF(4, UnitSquare), 2, 4, 3)) == [[(1, 0)], [(2, 0), (4, 3)], [(3, 0)]]
        # Rotations wrap around the number of sites.
        @test plain(Crafts.merge_vertices(VMF(2, UnitTriangle), 1, 2, 5)) == [[(1, 0), (2, 2)]]

        @test_throws ArgumentError Crafts.merge_vertices(vmf, 2, 1, 0)
        @test_throws ArgumentError Crafts.merge_vertices(vmf, 2, 2, 0)
    end

    @testset "isdescendable" begin
        # Against brute force: every form reachable from the initial one, and for each the forms
        # reachable from it by one or more merges.
        for species in (UnitSquare, UnitTriangle)
            n = nsites(species)
            initial = VMF(4, species)
            forms = Dict(plain(initial) => initial)
            children = Dict{Any,Vector{Any}}()
            frontier = [initial]
            while !isempty(frontier)
                next = typeof(initial)[]
                for f in frontier
                    m = length(f.vertices)
                    kids = children[plain(f)] = Any[]
                    for i in 1:m, j in i+1:m, Δψ in 0:n-1
                        c = Crafts.merge_vertices(f, i, j, Δψ)
                        push!(kids, plain(c))
                        haskey(forms, plain(c)) || (forms[plain(c)] = c; push!(next, c))
                    end
                end
                frontier = next
            end

            reachable = Dict{Any,Set{Any}}()
            function descendants(k)
                haskey(reachable, k) && return reachable[k]
                s = Set{Any}()
                for c in children[k]
                    push!(s, c)
                    union!(s, descendants(c))
                end
                return reachable[k] = s
            end

            @test all(Crafts.isdescendable(a, b) == (ka ∈ descendants(kb)) for (kb, b) in forms, (ka, a) in forms)
        end

        a = VMF([[(1, 0), (2, 1)]], 4)
        b = VMF(2, UnitSquare)
        @test Crafts.isdescendable(a, b)
        @test !Crafts.isdescendable(b, a)
        # Strict: a form is not its own descendant.
        @test !Crafts.isdescendable(b, b)
        # Same grouping, but species 2 sits at a different rotation relative to species 1.
        c = VMF([[(1, 0), (2, 1)], [(3, 0)]], 4)
        @test Crafts.isdescendable(VMF([[(1, 0), (2, 1), (3, 2)]], 4), c)
        @test !Crafts.isdescendable(VMF([[(1, 0), (2, 2), (3, 2)]], 4), c)
        # A different number of sites, or of initial species.
        @test !Crafts.isdescendable(a, VMF(2, UnitTriangle))
        @test !Crafts.isdescendable(a, VMF(3, UnitSquare))
    end

    ### Isomorphic rule sets for the maps: two generations of merges of the 2x2 block of four distinct
    ### squares from the bindingruledesign tests, built alongside their forms.
    n = nsites(UnitSquare)
    rules = BindingRules([1 2 2 4; 2 3 3 1; 3 4 4 2; 4 1 1 3], UnitSquare)
    systems = [(rules, VMF(4, UnitSquare))]
    generation = copy(systems)
    for _ in 1:2
        generation = [(Crafts.merge_species(r, i, j, Δψ), Crafts.merge_vertices(v, i, j, Δψ))
                      for (r, v) in generation for i in 1:nspecies(r) for j in i+1:nspecies(r) for Δψ in 0:n-1]
        append!(systems, generation)
    end
    graphs = [graphrep(r) for (r, _) in systems]
    perms = [Crafts.canonize!(g) for g in graphs]
    getmap(a, b) = Crafts.get_map_from_perm(systems[a][2], systems[b][2], perms[a], perms[b])
    identity_map(m) = [(vertex=k, rotation=0) for k in 1:m]

    # Bonds as (species, site, species, site), one color per site.
    function bondset(r)
        imat = interactionmatrix(r)
        return Set((divrem(c1 - 1, n) .+ 1..., divrem(c2 - 1, n) .+ 1...)
                   for c1 in axes(imat, 1), c2 in axes(imat, 2) if imat[c1, c2])
    end
    bondsets = [bondset(r) for (r, _) in systems]

    # Whether `mp`, taking vertices of the rule set with bonds `bonds_b` to those of the one with bonds
    # `bonds_a`, is a bijection on vertices carrying the first set of bonds onto the second. Site k' of x
    # is site k' - z of y.
    function maps_bonds(mp, bonds_a, bonds_b)
        sort([e.vertex for e in mp]) == 1:length(mp) || return false
        f(x, k) = (mp[x].vertex, mod1(k - mp[x].rotation, n))
        return Set((f(s1, k1)..., f(s2, k2)...) for (s1, k1, s2, k2) in bonds_b) == bonds_a
    end

    pairs = [(a, b) for a in eachindex(systems), b in eachindex(systems) if graphs[a] == graphs[b] && a != b]

    @testset "get_map_from_perm" begin
        @test !isempty(pairs)
        @test all(maps_bonds(getmap(a, b), bondsets[a], bondsets[b]) for (a, b) in pairs)

        @test getmap(1, 1) == identity_map(4)

        @test_throws ArgumentError getmap(1, 2)
        @test_throws ArgumentError Crafts.get_map_from_perm(systems[1][2], VMF(4, UnitTriangle), perms[1], perms[1])
    end

    @testset "invert_map and compose_maps" begin
        # Vertex 1 is vertex 2 at rotation 1, so vertex 2 is vertex 1 at rotation -1 ≡ 3, and so on.
        φ = [(vertex=2, rotation=1), (vertex=3, rotation=0), (vertex=1, rotation=2)]
        @test Crafts.invert_map(φ, 4) == [(vertex=3, rotation=2), (vertex=1, rotation=3), (vertex=2, rotation=0)]
        @test Crafts.invert_map(Crafts.invert_map(φ, 4), 4) == φ

        # Applying ψ and then φ: vertex 1 goes to 3 at rotation 3, and then to 1 at rotation 3 + 2.
        ψ = [(vertex=3, rotation=3), (vertex=1, rotation=0), (vertex=2, rotation=1)]
        @test Crafts.compose_maps(φ, ψ, 4) == [(vertex=1, rotation=1), (vertex=2, rotation=1), (vertex=3, rotation=1)]
        @test Crafts.compose_maps(φ, identity_map(3), 4) == φ
        @test Crafts.compose_maps(identity_map(3), φ, 4) == φ

        # A map composed with its inverse, in either order, is the identity.
        @test all(pairs) do (a, b)
            mp = getmap(a, b)
            inv = Crafts.invert_map(mp, n)
            m = length(mp)
            Crafts.compose_maps(mp, inv, n) == identity_map(m) && Crafts.compose_maps(inv, mp, n) == identity_map(m)
        end

        # Both agree with mapping the rule sets directly, which itself carries the bonds over.
        @test all(Crafts.invert_map(getmap(a, b), n) == getmap(b, a) for (a, b) in pairs)
        groups = [findall(==(g), graphs) for g in unique(graphs)]
        triples = [(a, b, c) for grp in groups for a in grp, b in grp, c in grp if a != b != c != a]
        @test !isempty(triples)
        @test all(triples) do (a, b, c)
            mp = Crafts.compose_maps(getmap(a, b), getmap(b, c), n)
            mp == getmap(a, c) && maps_bonds(mp, bondsets[a], bondsets[c])
        end
    end

    @testset "create_merge_isomorphs" begin
        # The rule set a form describes, found among the merges of a system's rules.
        function rules_of(vmf, a)
            r, v = systems[a]
            for i in 1:nspecies(r), j in i+1:nspecies(r), Δψ in 0:n-1
                Crafts.merge_vertices(v, i, j, Δψ).vertices == vmf.vertices && return Crafts.merge_species(r, i, j, Δψ)
            end
            error("form is not a merge of system $a")
        end

        # Every merge of every system, carried over to every other system isomorphic to it: the merges
        # must match those of the rule sets, with maps carrying the bonds of one onto the other.
        nchecked = 0
        allcorrect = true
        for b in eachindex(systems)
            others = [a for (a, b′) in pairs if b′ == b]
            isempty(others) && continue
            (r, v) = systems[b]
            for i in 1:nspecies(r), j in i+1:nspecies(r), Δψ in 0:n-1
                merged = bondset(Crafts.merge_species(r, i, j, Δψ))
                isomorphs, isomaps = Crafts.create_merge_isomorphs(v, i, j, Δψ, [systems[a][2] for a in others],
                                                                    [getmap(a, b) for a in others])
                for (a, iso, isomap) in zip(others, isomorphs, isomaps)
                    allcorrect &= maps_bonds(isomap, bondset(rules_of(iso, a)), merged)
                    nchecked += 1
                end
            end
        end
        @test nchecked > 0
        @test allcorrect

        # Against the form itself, under the identity, the merge is the same and the map stays the identity.
        vmf = VMF([[(1, 0), (5, 3)], [(2, 0), (3, 2)], [(4, 0)]], 4)
        isomorphs, isomaps = Crafts.create_merge_isomorphs(vmf, 1, 3, 1, [vmf], [identity_map(3)])
        @test only(isomorphs).vertices == Crafts.merge_vertices(vmf, 1, 3, 1).vertices
        @test only(isomaps) == identity_map(2)

        @test Crafts.create_merge_isomorphs(vmf, 1, 3, 1, VMF[], Vector{Crafts.VertexAtRotation}[]) == (VMF[], [])
        @test_throws ArgumentError Crafts.create_merge_isomorphs(vmf, 3, 1, 1, [vmf], [identity_map(3)])
        @test_throws ArgumentError Crafts.create_merge_isomorphs(vmf, 1, 3, 1, [vmf], typeof(identity_map(3))[])
    end
end
