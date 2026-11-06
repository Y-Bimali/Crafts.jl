### Drawing of a planar `BindingRules` of polygons as a diagram: each species drawn once, as one line per
### edge in the color of its site, with a curve joining every two binding sites the rules let bind. Several
### rule sets can be drawn side by side, each on its own terms.
###
### Needs CairoMakie and GLMakie in the active environment.

using Roly
using LinearAlgebra: Diagonal, Symmetric, diagind, dot, eigen, norm, pinv
import Random
using CairoMakie
import GLMakie

# As in plottargettree.jl: GLMakie shows the figures the REPL returns, CairoMakie only saves files.
GLMakie.activate!()

# Roly's Makie extension, whose palette is reused so that the sites have the colors `render` gives them. It
# is not public, so a change in Roly can break this file.
const _ROLY_MAKIE = Base.get_extension(Roly, :MakieExt)

# Geometry, in units of the largest radius among the species: the gap between neighbouring species, the
# reach of a self-loop out of its site and its width, the least distance a curve runs along the normal of
# each of its sites before it turns, the sideways bows a curve is tried with to get around the species in
# its way, and how far out of its site a site's color label sits.
const RULES_SPACING = 1.6
const RULES_LOOP_REACH = 0.9
const RULES_LOOP_WIDTH = 0.7
const RULES_MIN_HANDLE = 0.5
const RULES_BOWS = (0.0, 0.6, -0.6, 1.2, -1.2, 2.0, -2.0, 3.0, -3.0)
const RULES_LABEL_OFFSET = 0.2
const RULES_CURVE_SAMPLES = 48

# What a routed curve pays for each of its samples inside a species, for each curve it crosses, and for
# each unit of bow, in units of the largest species radius. Running through a species is the worst.
const RULES_COST_INSIDE = 1.0
const RULES_COST_CROSSING = 2.0
const RULES_COST_BOW = 0.3

# How much each unit of turning, in radians, a species needs costs against how badly the species face their
# partners, in choosing the way round the layout is drawn: half a turn costs about as much as a bonding site
# facing 35° away from its partner.
const RULES_TURN_COST = 0.06

# The height of a species' label, in units of the radius of the empty disk about the middle of the species:
# as large as leaves room in it, turned any way, for one underlined digit when no species has more, and for
# two, the first a 1, when some does. And, in units of that height, the width per digit of the underline that
# marks a label that reads as another upside down and how far below the label's middle it runs; and the
# underline's width in pixels.
const RULES_SPECIES_LABEL_SIZE_1DIGIT = 1.6
const RULES_SPECIES_LABEL_SIZE_2DIGITS = 1.3
const RULES_UNDERLINE_WIDTH = 0.5
const RULES_UNDERLINE_DROP = 0.5
const RULES_UNDERLINE_LINEWIDTH = 1.5

# The color of an edge that cannot bond: darker than Roly's, which is meant for a filled band and all but
# vanishes as a line on white.
const POLYGON_INERT_COLOR = RGBf(0.78, 0.78, 0.78)

# What each digit reads as upside down, for the digits that read as one.
const _UPSIDEDOWN_DIGITS = Dict('0' => '0', '1' => '1', '6' => '9', '8' => '8', '9' => '6')

# Whether `label`, a number, reads upside down as a different number, as 6, 9, 16 and 66 do.
function _upsidedown_ambiguous(label::AbstractString)
    all(in(keys(_UPSIDEDOWN_DIGITS)), label) || return false
    turned = String([_UPSIDEDOWN_DIGITS[c] for c in reverse(label)])
    return turned != label && !startswith(turned, '0')
end

# The radius of a disk about the origin of the polygon's frame that holds the polygon.
_speciesradius(spcs::PolygonParticleSpecies) = Float64(maximum(norm, spcs.corners))

### Layout of the species

# The number of bonds between each two species of `rules`, which bonds within a species do not count in.
function _speciesgraph(rules)
    n = nspecies(rules)
    A = zeros(Int, n, n)
    for (i, j) in Roly.bonded_species(rules)
        i == j && continue
        A[i, j] += 1
        A[j, i] += 1
    end
    return A
end

# The number of bonds on the shortest path between each two species of the graph `A`, infinite between
# species in separate components.
function _graphdistances(A)
    n = size(A, 1)
    D = [i == j ? 0.0 : A[i, j] > 0 ? 1.0 : Inf for i in 1:n, j in 1:n]
    for k in 1:n, i in 1:n, j in 1:n
        D[i, j] = min(D[i, j], D[i, k] + D[k, j])
    end
    return D
end

# The connected components of the graph whose distances are `D`, each as the sorted species in it, in the
# order of their first species.
_components(D) = unique([findall(isfinite, D[i, :]) for i in axes(D, 1)])

# The stress of the positions `X`, one per row, against the target distances `D`: how far their distances are
# from the targets, with each pair weighted by `1/D²` so that neighbours count the most.
_stress(X, D) = sum(((norm(X[i, :] - X[j, :]) - D[i, j]) / D[i, j])^2 for i in axes(X, 1), j in axes(X, 1) if i < j;
                    init=0.0)

# Positions, one per row, whose distances come as close to `D` as stress majorization gets them from `X`.
function _majorize(X, D; iterations=300)
    n = size(X, 1)
    W = [i == j ? 0.0 : 1 / D[i, j]^2 for i in 1:n, j in 1:n]
    Vinv = pinv(Diagonal(vec(sum(W; dims=2))) - W)
    for _ in 1:iterations
        B = map(CartesianIndices(W)) do ij
            i, j = Tuple(ij)
            dist = norm(X[i, :] - X[j, :])
            return i == j || dist == 0 ? 0.0 : -W[i, j] * D[i, j] / dist
        end
        B[diagind(B)] .= -vec(sum(B; dims=2))
        Xnew = Vinv * B * X
        maximum(abs, Xnew - X) < 1e-6 && return Xnew
        X = Xnew
    end
    return X
end

# The positions `X`, one per row, pushed apart pair by pair until no two are closer than 1, as happens where
# more species bond to each other than fit around one another in the plane.
function _separate(X; iterations=100)
    X = copy(X)
    for _ in 1:iterations
        moved = false
        for i in axes(X, 1), j in axes(X, 1)
            i < j || continue
            v = X[j, :] - X[i, :]
            dist = norm(v)
            dist >= 1 - 1e-9 && continue
            u = dist == 0 ? [1.0, 0.0] : v / dist
            X[i, :] -= (1 - dist) / 2 * u
            X[j, :] += (1 - dist) / 2 * u
            moved = true
        end
        moved || break
    end
    return X
end

# Positions, one per row, for the connected bond graph whose distances are `D`, a unit apart for any two
# species that bond: the best of a few runs of stress majorization, from a circle and from scattered starts,
# centred on the origin and turned so that it is widest along x.
function _componentlayout(D; nstarts=8)
    n = size(D, 1)
    n == 1 && return zeros(1, 2)
    rng = Random.Xoshiro(1)
    circle = [[cos(2π * i / n) for i in 1:n] [sin(2π * i / n) for i in 1:n]] .* (n / 2π)
    starts = [circle, (randn(rng, n, 2) .* (n / 2π) for _ in 2:nstarts)...]
    X = _separate(argmin(X -> _stress(X, D), (_majorize(X0, D) for X0 in starts)))
    X = X .- sum(X; dims=1) ./ n
    v = eigen(Symmetric(X' * X)).vectors[:, end]
    return X * [v[1] -v[2]; v[2] v[1]]
end

# Positions for the species with bond graph `A`, `d` apart for any two that bond. Each connected group of
# species is laid out on its own, as `_componentlayout` does, and the groups are set side by side, largest
# first, one spacing apart and centred on the x axis.
function _graphlayout(A, d)
    D = _graphdistances(A)
    comps = sort(_components(D); by=c -> -length(c))
    out = Vector{Tuple{Float64,Float64}}(undef, size(A, 1))
    x = 0.0
    for comp in comps
        X = _componentlayout(D[comp, comp])
        lo, hi = minimum(X[:, 1]), maximum(X[:, 1])
        for (k, i) in enumerate(comp)
            out[i] = (d * (x + X[k, 1] - lo), d * X[k, 2])
        end
        x += hi - lo + 1
    end
    return out
end

# The poses of the species: given as a vector of poses or of positions, or one of the layouts `:graph`,
# `:row`, `:grid` (`layoutcols` columns) and `:circle`. `:graph` sets the species apart by the bonds between
# them, as `_graphlayout` does with the bond graph `A`; `:auto` is the same.
function _speciesposes(n::Integer, layout, d::Real; layoutcols=nothing, A=nothing)
    layout isa AbstractVector && return [p isa Pose ? p : Pose(SVector{2,Float64}(p...), Angle2d(0.0))
                                         for p in layout]
    layout == :auto && (layout = :graph)
    xs = if layout == :graph
        _graphlayout(A, d)
    elseif layout == :row
        [(d * (i - 1), 0.0) for i in 1:n]
    elseif layout == :grid
        nc = something(layoutcols, ceil(Int, sqrt(n)))
        [(d * mod(i - 1, nc), -d * div(i - 1, nc)) for i in 1:n]
    elseif layout == :circle
        R = n == 1 ? 0.0 : d / (2sin(π / n))
        [(R * cos(π / 2 - 2π * (i - 1) / n), R * sin(π / 2 - 2π * (i - 1) / n)) for i in 1:n]
    else
        throw(ArgumentError("unknown layout $(repr(layout)), use :auto, :graph, :row, :grid, :circle or a vector"))
    end
    return [Pose(SVector{2,Float64}(x...), Angle2d(0.0)) for x in xs]
end

### Orientation of the species

# The outward unit normal of site `i` of `spcs` in the species' frame.
function _sitenormal(spcs, i)
    n = bindingsite(spcs, i).pose.psi * SVector(1.0, 0.0)
    return n / norm(n)
end

# The position and the outward unit normal, both in the plane of the axis, of site `i` of `spcs` at `pose`.
function _sitepoint(spcs, pose, i)
    sp = pose * bindingsite(spcs, i).pose
    return Point2f(sp.x...), Vec2f((pose.psi * _sitenormal(spcs, i))...)
end

# How badly the bonding sites of species `i` face the species their partners belong to when it is turned by
# `θ`: for each such site, one minus the cosine of the angle between its normal and the direction of the
# partner's species. Zero if the species bonds only to itself.
function _facingcost(spcs, poses, curves, i, θ)
    ends = [(s.site, o.species) for (a, b, _) in curves for (s, o) in ((a, b), (b, a))
            if s.species == i && o.species != i]
    return sum(ends; init=0.0) do (site, j)
        d = poses[j].x - poses[i].x
        return 1 - dot(Angle2d(θ) * poses[i].psi * _sitenormal(spcs[i], site), d / norm(d))
    end
end

# The angle, a multiple of 15°, by which each species is best turned so that its bonding sites face the
# species their partners belong to, and the total of `_facingcost` it leaves. The middles of the species do
# not move, so each species is turned on its own. Of equally good angles the smallest wins.
function _autorotations(spcs, poses, curves)
    angles = 2π .* (0:23) ./ 24
    rotations = map(eachindex(spcs)) do i
        argmin(θ -> _facingcost(spcs, poses, curves, i, θ) + 1e-3 * min(θ, 2π - θ), angles)
    end
    return rotations, sum(i -> _facingcost(spcs, poses, curves, i, rotations[i]), eachindex(spcs))
end

# The layout `poses` the way round, among itself turned by half a turn and mirrored about either axis
# through its middle, which all keep it as wide, that lets the species face their partners best and, of the
# ways that do so about as well, turn the least: a graph layout fixes only where species sit relative to
# each other, not which way round the whole drawing goes. Returns the poses and the turns of the species.
function _orientlayout(spcs, poses, curves)
    c = sum(p -> p.x, poses) / length(poses)
    best = nothing
    for (sx, sy) in ((1, 1), (-1, -1), (-1, 1), (1, -1))
        cand = [Pose(c + SVector(sx * (p.x - c)[1], sy * (p.x - c)[2]), p.psi) for p in poses]
        rotations, facing = _autorotations(spcs, cand, curves)
        cost = facing + RULES_TURN_COST * sum(θ -> min(θ, 2π - θ), rotations)
        (isnothing(best) || cost < best[1] - 1e-9) && (best = (cost, cand, rotations))
    end
    return best[2], best[3]
end

### Curves between the sites

# Points along the cubic Bézier curve with end points `p0`, `p3` and control points `p1`, `p2`.
function _bezier(p0, p1, p2, p3; n=RULES_CURVE_SAMPLES)
    return map(range(0, 1; length=n)) do t
        s = 1 - t
        s^3 * p0 + 3s^2 * t * p1 + 3s * t^2 * p2 + t^3 * p3
    end
end

# The loop out of and back into the site at `p` with normal `np`.
function _selfloop(p, np, r)
    t = Vec2f(-np[2], np[1])
    return _bezier(p, p + RULES_LOOP_REACH * r * (np + RULES_LOOP_WIDTH * t),
                   p + RULES_LOOP_REACH * r * (np - RULES_LOOP_WIDTH * t), p)
end

# The curve from the site at `p` with normal `np` to the one at `q` with normal `nq`, leaving and entering
# each along its normal, bowed sideways by `bow`.
function _bondcurve(p, np, q, nq, r, bow)
    h = max(RULES_MIN_HANDLE * r, 0.4norm(q - p))
    u = (q - p) / norm(q - p)
    side = bow * r * Vec2f(-u[2], u[1])
    return _bezier(p, p + h * np + side, q + h * nq + side, q)
end

# Whether the segments `a1`–`a2` and `b1`–`b2` cross.
function _crosses(a1, a2, b1, b2)
    c(o, p, q) = (p[1] - o[1]) * (q[2] - o[2]) - (p[2] - o[2]) * (q[1] - o[1])
    return c(a1, a2, b1) * c(a1, a2, b2) < 0 && c(b1, b2, a1) * c(b1, b2, a2) < 0
end

_ncrossings(pts, other) =
    count(_crosses(pts[i], pts[i + 1], other[j], other[j + 1]) for i in 1:(length(pts) - 1), j in 1:(length(other) - 1))

# What `pts` costs as a route: its samples inside the disks holding the species, given as `(middle, radius)`,
# the curves already routed that it crosses, and its bow. Samples within `reach` of either end of the curve
# are free, since a curve has to cross the disk of a species to get out of one of its sites.
function _routecost(pts, bow, disks, routed, reach)
    p, q = first(pts), last(pts)
    inside = count(norm(x - p) >= reach && norm(x - q) >= reach && any(norm(x - c) < ρ for (c, ρ) in disks)
                   for x in pts)
    crossings = isempty(routed) ? 0 : sum(other -> _ncrossings(pts, other), routed)
    return RULES_COST_INSIDE * inside + RULES_COST_CROSSING * crossings + RULES_COST_BOW * abs(bow)
end

# The pairs of sites joined by a curve, as `(a, b, k)` with `a` and `b` the `SpeciesSiteLoc`s and `k` the
# index of their pair of colors in `bonded_colors`. Unless `allsites`, each pair of colors gets one curve,
# between the first site of each color: sites that share a color are the same as far as the rules go.
function _bondcurves(rules; allsites=false)
    out = Tuple{SpeciesSiteLoc,SpeciesSiteLoc,Int}[]
    for (k, (c1, c2)) in enumerate(Roly.bonded_colors(rules))
        s1, s2 = Roly.sitesofcolor(rules, c1), Roly.sitesofcolor(rules, c2)
        if !allsites
            push!(out, (first(s1), first(s2), k))
            continue
        end
        for a in s1, b in s2
            c1 == c2 && b < a && continue
            push!(out, (a, b, k))
        end
    end
    return out
end

# The curves of `curves`, routed between the species at `poses` as `plot_polygonrules!` describes: the
# self-loops first, since they cannot move, and then the shortest curves, which have the fewest ways around
# what is in their way.
function _routecurves(spcs, poses, curves, r; route)
    disks = [(Point2f(p.x...), Float32(_speciesradius(s))) for (s, p) in zip(spcs, poses)]
    ends = map(curves) do (a, b, _)
        (_sitepoint(spcs[a.species], poses[a.species], a.site)..., _sitepoint(spcs[b.species], poses[b.species], b.site)...)
    end
    order = sortperm(ends; by=e -> e[1] ≈ e[3] ? -1.0 : Float64(norm(e[3] - e[1])))
    routed = Vector{Vector{Point2f}}(undef, length(curves))
    for (m, c) in enumerate(order)
        p, np, q, nq = ends[c]
        routed[c] = if p ≈ q
            _selfloop(p, np, r)
        elseif route
            others = routed[order[1:(m - 1)]]
            bow = argmin(b -> _routecost(_bondcurve(p, np, q, nq, r, b), b, disks, others, RULES_MIN_HANDLE * r), RULES_BOWS)
            _bondcurve(p, np, q, nq, r, bow)
        else
            _bondcurve(p, np, q, nq, r, 0.0)
        end
    end
    return routed
end

### Drawing

# Draw `spcs` at `pose` as one line per edge, edge `j` from corner `j` to corner `j + 1`, which is where Roly
# puts site `j`, in `colors[j]`. Each line stops short of the corners by `cornergap` of the edge, so that the
# sites read as separate.
function _drawpolygon!(ax, spcs, pose, colors; edgewidth, cornergap)
    corners = [Point2f((pose * SVector{2,Float64}(c...))...) for c in spcs.corners]
    n = length(corners)
    for j in 1:n
        a, b = corners[j], corners[mod1(j + 1, n)]
        lines!(ax, [a + cornergap * (b - a), b - cornergap * (b - a)]; color=colors[j], linewidth=edgewidth)
    end
end

# Write `i` in the middle of the species at `pose`, `h` high and turned with it, underlined if it reads as
# another number upside down.
function _drawspecieslabel!(ax, i, pose, h)
    text!(ax, Point2f(pose.x...); text=string(i), align=(:center, :center), markerspace=:data, fontsize=h,
          rotation=rotation_angle(pose.psi))
    if _upsidedown_ambiguous(string(i))
        w = RULES_UNDERLINE_WIDTH * h * ndigits(i)
        ends = [Point2f((pose * SVector(x, -RULES_UNDERLINE_DROP * h))...) for x in (-w / 2, w / 2)]
        lines!(ax, ends; color=:black, linewidth=RULES_UNDERLINE_LINEWIDTH)
    end
end

"""
    plot_polygonrules!(ax, rules::BindingRules; layout=:auto, layoutcols=nothing, rotations=:auto,
                       spacing=RULES_SPACING, allsites=false, route=true, bondcolor=:gradient, linewidth=3,
                       edgewidth=5, cornergap=0.06, inertcolor=POLYGON_INERT_COLOR, specieslabels=true,
                       colorlabels=false, fontsize=14)

Draw planar `rules`, whose species must all be polygons, onto `ax`: every species once, as one line per
edge in the color Roly gives the edge's site, and a curve joining every two binding sites the rules let
bind. A site that binds a copy of itself gets a loop.

`layout` places the species: `:graph` (the default, also `:auto`), `:row`, `:grid` (in `layoutcols`
columns), `:circle` or a vector of one `Pose` or position per species. `:graph` lays out the graph whose
nodes are the species, joined where any of their sites bond, so that the species are about as far apart as
the fewest bonds between them: species that bond lie next to each other, and separate groups apart.
`spacing` is the gap between neighbouring species, in units of the largest species radius.

`rotations` turns each species on top of the orientation the layout gives it: by a given angle each, in
radians, by none if `nothing`, or, with `:auto`, so that its bonding sites face their partners' species.
With `:auto` and the graph layout, the layout is also turned by half a turn or mirrored, if that lets the
species face their partners better or turn less.

Unless `allsites`, sites that share a color are drawn joined only once, through the first of them.

With `route`, each curve is bowed sideways as far as it takes to keep it out of the species and to cross
as few curves as it can; without it every curve takes the direct route. `linewidth` is the curves' width in
pixels.

`bondcolor` is `:gradient`, a curve shading from the color of one site to that of the other, `:bond`, a
color of its own for each pair of colors in the rules, or any color for all curves.

The edges are `edgewidth` pixels wide, those whose site cannot bond drawn in `inertcolor`. Each stops short
of the corners by `cornergap` of its length, so that the sites read as separate.

`specieslabels` writes the index of each species in its middle, turned with the species, so that it shows
which way round the species is drawn and hence which site is which. A number that reads as another upside
down, such as 6 or 9, is underlined. `colorlabels` writes just outside each site its color in the rules,
`fontsize` pixels high.
"""
function plot_polygonrules!(ax, rules::BindingRules; layout=:auto, layoutcols=nothing, rotations=:auto,
                            spacing=RULES_SPACING, allsites=false, route=true, bondcolor=:gradient,
                            linewidth=3, edgewidth=5, cornergap=0.06, inertcolor=POLYGON_INERT_COLOR,
                            specieslabels=true, colorlabels=false, fontsize=14)
    dimension(rules) == 2 || throw(ArgumentError("only planar rules can be drawn"))
    spcs = species(rules)
    all(s -> s isa PolygonParticleSpecies, spcs) ||
        throw(ArgumentError("plot_polygonrules! draws rules of polygons only"))

    n = length(spcs)
    r = maximum(_speciesradius, spcs)
    poses = _speciesposes(n, layout, (2 + spacing) * r; layoutcols, A=_speciesgraph(rules))
    length(poses) == n || throw(ArgumentError("layout gives $(length(poses)) poses for $n species"))
    curves = _bondcurves(rules; allsites)
    if rotations == :auto
        poses, rotations = layout in (:auto, :graph) ? _orientlayout(spcs, poses, curves) :
                                                       (poses, first(_autorotations(spcs, poses, curves)))
    end
    isnothing(rotations) || (poses = [Pose(p.x, Angle2d(θ) * p.psi) for (p, θ) in zip(poses, rotations)])

    sitecolors = [_ROLY_MAKIE._resolve_colors(spcs[i], i, rules, nothing)[3] for i in 1:n]
    for i in 1:n
        colors = [isinert(rules, SpeciesSiteLoc(i, j)) ? inertcolor : sitecolors[i][j] for j in 1:nsites(spcs[i])]
        _drawpolygon!(ax, spcs[i], poses[i], colors; edgewidth, cornergap)
    end

    routed = _routecurves(spcs, poses, curves, r; route)
    bondpalette = Makie.wong_colors()
    for (c, (a, b, k)) in enumerate(curves)
        pts = routed[c]
        curvecolor = if bondcolor == :gradient
            ca, cb = sitecolors[a.species][a.site], sitecolors[b.species][b.site]
            [Makie.Colors.weighted_color_mean(1 - t, RGBf(ca), RGBf(cb)) for t in range(0, 1; length=length(pts))]
        elseif bondcolor == :bond
            bondpalette[mod1(k, length(bondpalette))]
        else
            bondcolor
        end
        lines!(ax, pts; color=curvecolor, linewidth, linecap=:round, joinstyle=:round)
        scatter!(ax, [first(pts), last(pts)]; color=:black, markersize=linewidth + 4)
    end

    if specieslabels
        # Sized to the particle rather than the screen, so that it stays inside the particle however far the
        # figure is zoomed, in the room the edges leave, which are lines a few pixels wide, a little short of
        # the inradius. One size for all species with the same room, so that 9 and 10 match.
        labelsize = n <= 9 ? RULES_SPECIES_LABEL_SIZE_1DIGIT : RULES_SPECIES_LABEL_SIZE_2DIGITS
        for (i, p) in enumerate(poses)
            _drawspecieslabel!(ax, i, p, labelsize * 0.9spcs[i].rmin)
        end
    end
    if colorlabels
        # Beside the curve leaving the site rather than on it.
        for i in 1:n, j in 1:nsites(spcs[i])
            p, np = _sitepoint(spcs[i], poses[i], j)
            text!(ax, p + RULES_LABEL_OFFSET * r * (np + Vec2f(-np[2], np[1]));
                  text=string(color(rules, SpeciesSiteLoc(i, j))), align=(:center, :center), fontsize=0.85fontsize)
        end
    end
    return ax
end

### Figures

# The `(row, column)` of the plot of each of `m` rule sets, in order, and the number of rows and columns of
# the grid they lie in: the cells of `arrangement` that are not zero, row by row, if it is given, and
# otherwise the cells of `ncols` columns, row by row.
function _plotcells(m::Integer, ncols::Integer, arrangement)
    if isnothing(arrangement)
        ncols >= 1 || throw(ArgumentError("ncols must be at least 1, not $ncols"))
        return [(div(k - 1, ncols) + 1, mod(k - 1, ncols) + 1) for k in 1:m], (cld(m, ncols), min(m, ncols))
    end
    cells = [(i, j) for i in axes(arrangement, 1) for j in axes(arrangement, 2) if !iszero(arrangement[i, j])]
    length(cells) == m ||
        throw(ArgumentError("arrangement has $(length(cells)) cells to plot in, for $m rule sets"))
    return cells, size(arrangement)
end

# The room, in units of length, left around each plot's drawing, so that the curves' width is not cut off at
# its edges.
const RULES_PLOT_MARGIN = 0.2

"""
    plot_polygonrules(rulesets::AbstractVector{<:BindingRules}; titles=nothing, titlesize=16, ncols=3,
                      arrangement=nothing, unitsize=60, width=nothing, filename=nothing, kwargs...)

Draw each rule set of `rulesets` as [`plot_polygonrules!`](@ref) does with `kwargs`, in a grid of plots of
one figure. Each plot is drawn on its own terms, so that the same species numbers and colors come up in
several of them, but all to one scale, so that a unit square is as large in each.

The plots fill the grid in order, row by row: `ncols` columns of it, or, if `arrangement` is given, the
cells of the matrix `arrangement` that are true or not zero, so that `[0 1 0; 1 1 1; 0 1 0]` sets five plots
out in a cross, and an empty row or column of it leaves a gap.

`titles` holds the title of each plot, in the order of `rulesets`; a title of `nothing` leaves its plot
without one, and so does `titles` of `nothing` all. `titlesize` is their height in pixels. A title may be a
`LaTeXString`, `L"..."`, to be set in TeX.

The scale is `unitsize` pixels to a unit of length, the side of a unit square, unless `width` is given, in
which case it is what makes the plots `width` pixels wide all together, without the gaps between them. Each
column is as wide as its widest plot, each row as high as its highest; an empty one, as wide or high as the
others are on average.

Returns the figure, after saving it to `filename` with CairoMakie if given. The REPL shows it in a GLMakie
window unless the call ends in `;`, elsewhere `display` it to get one.
"""
function plot_polygonrules(rulesets::AbstractVector{<:BindingRules}; titles=nothing, titlesize=16,
                           ncols::Integer=3, arrangement=nothing, unitsize=60, width=nothing, filename=nothing,
                           kwargs...)
    m = length(rulesets)
    m >= 1 || throw(ArgumentError("no rule sets to plot"))
    isnothing(titles) || length(titles) == m ||
        throw(ArgumentError("$(length(titles)) titles for $m rule sets"))
    cells, (nrows, ncolumns) = _plotcells(m, ncols, arrangement)

    fig = Figure()
    axes_ = Vector{Axis}(undef, m)
    boxes = Vector{Rect2f}(undef, m)
    for (k, (i, j)) in enumerate(cells)
        title = isnothing(titles) ? "" : something(titles[k], "")
        ax = axes_[k] = Axis(fig[i, j]; aspect=DataAspect(), title, titlesize)
        hidedecorations!(ax)
        hidespines!(ax)
        plot_polygonrules!(ax, rulesets[k]; kwargs...)
        box = Makie.data_limits(ax.scene)
        boxes[k] = Rect2f(Vec2f(minimum(box)[1:2]) .- RULES_PLOT_MARGIN, Vec2f(widths(box)[1:2]) .+ 2RULES_PLOT_MARGIN)
    end

    # The widths of the columns and the heights of the rows, in units of length.
    extent(d, n, on) = map(1:n) do l
        maximum(widths(boxes[k])[d] for (k, cell) in enumerate(cells) if cell[on] == l; init=-Inf)
    end
    average(xs) = (filled = filter(isfinite, xs); replace(xs, -Inf => sum(filled) / length(filled)))
    colwidths = average(extent(1, ncolumns, 2))
    rowheights = average(extent(2, nrows, 1))
    isnothing(width) || (unitsize = width / sum(colwidths))

    # Each plot centred in its cell, with the cell's extent, so that every plot has the same scale.
    for (k, (i, j)) in enumerate(cells)
        c = minimum(boxes[k]) + widths(boxes[k]) / 2
        w, h = colwidths[j], rowheights[i]
        limits!(axes_[k], c[1] - w / 2, c[1] + w / 2, c[2] - h / 2, c[2] + h / 2)
    end
    for j in 1:ncolumns
        colsize!(fig.layout, j, Fixed(round(Int, unitsize * colwidths[j])))
    end
    for i in 1:nrows
        rowsize!(fig.layout, i, Fixed(round(Int, unitsize * rowheights[i])))
    end
    resize_to_layout!(fig)

    isnothing(filename) || save(filename, fig; backend=CairoMakie)
    return fig
end

"""
    plot_polygonrules(rules::BindingRules; title=string(rules), width=600, filename=nothing, kwargs...)

Draw `rules` as [`plot_polygonrules!`](@ref) does with `kwargs`, in a figure of its own, `width` pixels
wide and as high as the drawing needs, under `title`, or none if `nothing`.

Returns the figure, after saving it to `filename` with CairoMakie if given, as the method for several rule
sets does.
"""
plot_polygonrules(rules::BindingRules; title=string(rules), width=600, kwargs...) =
    plot_polygonrules([rules]; titles=[title], ncols=1, width, kwargs...)
