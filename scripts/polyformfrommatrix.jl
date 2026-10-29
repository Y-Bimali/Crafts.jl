using Roly
using Roly: ParticleSiteLoc, SpeciesSiteLoc, bindingsite, raise!

"""
    polyform_from_matrix(shape::AbstractMatrix{<:Integer}, species)

Return the `Polyform` that fills the cells of `shape` holding a one, built from `UnitSquare`s by
`_square_polyform_from_matrix` or from `UnitTriangle`s by `_triangle_polyform_from_matrix`,
whichever `species` is.

Throws an `ArgumentError` if `shape` holds anything but zeros and ones or holds no one, or if `species`
is neither `UnitSquare` nor `UnitTriangle`.
"""
function polyform_from_matrix(shape::AbstractMatrix{<:Integer}, species)
    all(x -> x == 0 || x == 1, shape) || throw(ArgumentError("the shape must hold only zeros and ones"))
    any(isone, shape) || throw(ArgumentError("the shape must hold at least one one"))

    species == UnitSquare && return _square_polyform_from_matrix(shape)
    species == UnitTriangle && return _triangle_polyform_from_matrix(shape)
    throw(ArgumentError("the species must be UnitSquare or UnitTriangle"))
end

"""
    _square_polyform_from_matrix(shape::AbstractMatrix{<:Integer})

Return the `Polyform` of `UnitSquare`s that fills the cells of `shape` holding a one, read as a picture:
row 1 is the top, column 1 the left.

The rules are a single species whose opposite sites bind. These are the only bonds squares of one
orientation can make, so every pair of adjacent squares is bonded and all of them share that orientation.
`shape` is taken to be checked already, as [`polyform_from_matrix`](@ref) does. Throws an
`ArgumentError` if its ones are not connected.
"""
function _square_polyform_from_matrix(shape::AbstractMatrix{<:Integer})
    cells = findall(isone, shape)
    rules = BindingRules([1 1 1 3;
                          1 2 1 4], UnitSquare)

    # UnitSquare sites face down, left, up and right. Each step is a neighbor's offset in the matrix,
    # the site of the placed square facing it, and the site of the new square facing back.
    steps = ((CartesianIndex(1, 0), 1, 3), (CartesianIndex(-1, 0), 3, 1),
             (CartesianIndex(0, 1), 4, 2), (CartesianIndex(0, -1), 2, 4))

    # Breadth first from the first square, so every new square attaches to one already placed.
    # Particles keep the order they are attached in, which is what `index` records.
    p = Roly.Polyform(rules, 1)
    index = Dict(first(cells) => 1)
    queue = [first(cells)]
    head = 1
    while head <= length(queue)
        cell = queue[head]
        head += 1
        for (offset, site, mate) in steps
            next = cell + offset
            checkbounds(Bool, shape, next) && isone(shape[next]) && !haskey(index, next) || continue
            attached = raise!(p, bindingsite(p, ParticleSiteLoc(index[cell], site)), SpeciesSiteLoc(1, mate))
            ismissing(attached) && error("could not attach the square at $(Tuple(next))")
            index[next] = nparticles(p)
            push!(queue, next)
        end
    end

    length(index) == length(cells) || throw(ArgumentError("the ones of the shape are not connected"))
    return p
end

"""
    _triangle_polyform_from_matrix(shape::AbstractMatrix{<:Integer})

Return the `Polyform` of `UnitTriangle`s that fills the cells of `shape` holding a one, read as a
picture on the triangular lattice: row 1 is the top, column 1 the left, and the cell in row `r` and
column `c` holds a triangle pointing up when `r + c` is even and pointing down otherwise. So
`[1 1 1; 1 1 1]` is a hexagon and `[1 0; 1 1; 0 1]` a parallelogram.

An up triangle shares its bottom edge with the down triangle below it, and each triangle shares its
slanted edges with its neighbors in the same row. Up triangles keep the species' own orientation and
down triangles are that rotated by half a turn. The binding rules are each edge to itself, meaning
every pair of adjacent triangles is bonded.
`shape` is taken to be checked already, as [`polyform_from_matrix`](@ref) does. Throws an
`ArgumentError` if its ones are not connected.
"""
function _triangle_polyform_from_matrix(shape::AbstractMatrix{<:Integer})
    cells = findall(isone, shape)
    rules = BindingRules(reduce(vcat, [1 i 1 i] for i in 1:3), UnitTriangle)

    # An up UnitTriangle has sites facing down, upper left and upper right; turning it by half a turn
    # makes them face up, lower right and lower left. Each step is a neighbor's offset in the matrix,
    # the site of the placed triangle facing it, and the site of the new triangle facing back. Joining
    # a site to its namesake is what puts the new triangle in the other of the two orientations.
    isup(cell) = iseven(sum(Tuple(cell)))
    upsteps = ((CartesianIndex(1, 0), 1, 1), (CartesianIndex(0, -1), 2, 2), (CartesianIndex(0, 1), 3, 3))
    downsteps = ((CartesianIndex(-1, 0), 1, 1), (CartesianIndex(0, -1), 3, 3), (CartesianIndex(0, 1), 2, 2))

    # Breadth first from an up triangle, which starts in the species' own orientation, so every new
    # triangle attaches to one already placed. Particles keep the order they are attached in, which is
    # what `index` records. Only a lone down triangle has no up triangle to start from.
    start = something(findfirst(isup, cells), 1)
    p = Roly.Polyform(rules, 1)
    index = Dict(cells[start] => 1)
    queue = [cells[start]]
    head = 1
    while head <= length(queue)
        cell = queue[head]
        head += 1
        for (offset, site, mate) in (isup(cell) ? upsteps : downsteps)
            next = cell + offset
            checkbounds(Bool, shape, next) && isone(shape[next]) && !haskey(index, next) || continue
            attached = raise!(p, bindingsite(p, ParticleSiteLoc(index[cell], site)), SpeciesSiteLoc(1, mate))
            ismissing(attached) && error("could not attach the triangle at $(Tuple(next))")
            index[next] = nparticles(p)
            push!(queue, next)
        end
    end

    length(index) == length(cells) || throw(ArgumentError("the ones of the shape are not connected"))
    return p
end
