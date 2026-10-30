### Plotting of the results of the design searches, kept apart from the searches themselves.
###
### Needs CairoMakie and GLMakie in the active environment.

using Crafts, Roly, Graphs
using CairoMakie
import GLMakie
using Base.Threads

# Whichever backend loaded last is active, and the REPL shows returned figures with it. Under CairoMakie
# that opens a PNG in external viewers, so GLMakie is made active here and CairoMakie only saves files.
GLMakie.activate!()

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

# Side, in the target's units, of the region shown by the target inset, the margin in those units below a
# planar target, the side of the inset in pixels, the gap in pixels between the inset and the line below it,
# and the width in pixels of the particles' outlines.
# TODO: make the size a parameter of `plot_target_tree`.
const TARGET_INSET_EXTENT = 5.5
const TARGET_INSET_MARGIN = 0.1
const TARGET_INSET_PX = 120
const TARGET_INSET_GAP = 4
const TARGET_INSET_STROKEWIDTH = 1

# Draw `target` as `render` would, but with thinner outlines, in a transparent inset whose bottom left corner
# sits just above the point (0, `y`) of `ax`, following it as `ax` moves or zooms. A planar target is centered
# horizontally on the middle of its bounding box and rests on the inset's bottom edge, one in space centered
# on the middle of its bounding box.
function _plot_target_inset!(ax::Axis, y, target)
    is3d = dimension(target) == 3
    bbox = lift(ax.scene.viewport, ax.finallimits) do viewport, lims
        bottom = viewport.origin[2] + viewport.widths[2] * (y - lims.origin[2]) / lims.widths[2]
        Rect2f(viewport.origin[1], bottom + TARGET_INSET_GAP, TARGET_INSET_PX, TARGET_INSET_PX)
    end
    placement = (; bbox, backgroundcolor=:transparent)
    inset = is3d ? Axis3(ax.parent; aspect=:data, placement...) : Axis(ax.parent; aspect=DataAspect(), placement...)
    hidedecorations!(inset)
    hidespines!(inset)
    # The polyformplot recipe passes no attributes on to the particles, so the extension's own drawing
    # function is called to set the width of their outlines.
    Base.get_extension(Roly, :MakieExt).plot_polyform!(inset, target; strokewidth=TARGET_INSET_STROKEWIDTH)

    box = Makie.data_limits(inset.scene)
    center = minimum(box) .+ widths(box) ./ 2
    lims = [(c - TARGET_INSET_EXTENT/2, c + TARGET_INSET_EXTENT/2) for c in center[1:(is3d ? 3 : 2)]]
    if !is3d
        bottom = minimum(box)[2] - TARGET_INSET_MARGIN
        lims[2] = (bottom, bottom + TARGET_INSET_EXTENT)
    end
    limits!(inset, Iterators.flatten(lims)...)
    return inset
end

# Gaps, in pixels, between the level axes and the n and count columns beside them, and between the title
# and the column headers.
const LABEL_COLGAP = 8
const TITLE_ROWGAP = 6

# The number of structures p and the number of species n of each rule set of the tree, the values its
# point is colored by under `colorby` with their label, and the tree's number of levels.
function _target_tree_data(rules::AbstractVector{<:BindingRules}, tree::AbstractGraph, colorby::Symbol)
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

    multiplicity = Dict{Tuple{Int,Int},Int}()
    foreach(np -> multiplicity[np] = get(multiplicity, np, 0) + 1, zip(ns, ps))
    colorvalues, colorlabel = if colorby == :bonds
        nbonds.(rules), L"bonds"
    elseif colorby == :maxsize
        maxsizes, "max. size"
    else
        [multiplicity[np] for np in zip(ns, ps)], "multiplicity"
    end

    nmin, nmax = extrema(ns)
    return (; ps, ns, colorvalues, colorlabel, nlevels=nmax - nmin + 1)
end

# The right end of the level axes: `maxpolys`, by default the largest p in `ps`.
_target_tree_xmax(maxpolys, ps) = something(maxpolys, max(maximum(ps), 1))

# Default height in pixels of a plot of a tree with `nlevels` levels.
_target_tree_height(nlevels) = 160 + 50nlevels

# Colormap and color range giving one color per integer from `vmin` to `vmax`, each centered on its value in
# a colorbar.
_integer_colors(colormap, vmin, vmax) =
    (; colormap=cgrad(colormap, vmax - vmin + 1; categorical=true), colorrange=(vmin - 0.5, vmax + 0.5))

# Colorbar at `pos` for the colors of `_integer_colors`.
function _integer_colorbar(pos, colormap, vmin, vmax, label)
    (; colormap, colorrange) = _integer_colors(colormap, vmin, vmax)
    return Colorbar(pos; label, colormap, limits=colorrange, ticks=_integer_ticks(vmin, vmax; maxticks=10))
end

# Draw the tree with the `_target_tree_data` `data` into `gl`, as described for `plot_target_tree`, its
# points colored with `colormap` over the integers from `vmin` to `vmax`. Without `title`, the column headers
# are the top row; without `showx`, the x axis has neither label nor tick labels; with `colorbar`, the
# colorbar is drawn in a fourth column. Returns the axis of the levels.
function _plot_target_tree!(gl::GridLayout, data, tree::AbstractGraph;
                            xmax, vmin, vmax, title, colormap, target, showx=true, colorbar=true)
    (; ps, ns, colorvalues, colorlabel, nlevels) = data
    nclipped = count(>(xmax), ps)
    nclipped > 0 && @warn "$nclipped rule sets allow more than maxpolys = $xmax structures and lie off the plot"

    # Level axes sit at y = 1, 2, ..., the root's on top, so that y = n - nmin + 1.
    nmin, nmax = extrema(ns)
    ys = ns .- nmin .+ 1

    header = isnothing(title) ? 1 : 2
    isnothing(title) || Label(gl[1, 1:(colorbar ? 4 : 3)], title; fontsize=20, font=:bold)
    Label(gl[header, 1], L"n_{spc}"; fontsize=18)
    Label(gl[header, 3], "count")

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
    showx || hidexdecorations!(ax; ticks=false, grid=false)

    for (row, n) in zip(levelrows, nmax:-1:nmin)
        Label(gl[row, 1], string(n))
        Label(gl[row, 3], string(count(==(n), ns)))
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

    linesegments!(ax, [Point2f(0, y) => Point2f(xmax, y) for y in 1:nlevels]; color=:gray50, linewidth=1)
    linesegments!(ax, [Point2f(ps[src(e)], ys[src(e)]) => Point2f(ps[dst(e)], ys[dst(e)]) for e in edges(tree)];
                  color=:black, linewidth=0.75)
    scatter!(ax, ps, ys; color=colorvalues, markersize=10, _integer_colors(colormap, vmin, vmax)...)
    colorbar && _integer_colorbar(gl[levelrows, 4], colormap, vmin, vmax, colorlabel)

    isnothing(target) || _plot_target_inset!(ax, nlevels, target)
    return ax
end

"""
    plot_target_tree(rules, tree; colorby=:bonds, maxpolys=nothing, title="Designable Search Tree",
                     filename=nothing, colormap=:viridis, figsize=nothing, target=nothing)

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

`target`, a `Polyform`, is drawn as by `render(target)`, with outlines 1 pixel wide, in an inset resting
on the left end of the root's axis, reaching up into the title. The inset is always in the same place and
shows a fixed region around the target, wide enough for a 4×4 structure of unit squares, so targets keep
their relative sizes across plots.

A `title` of `nothing` leaves the plot without one.

Returns the figure, after saving it to `filename` with CairoMakie if given. The REPL shows it in a
GLMakie window unless the call ends in `;`, elsewhere `display` it to get one.
"""
function plot_target_tree(rules::AbstractVector{<:BindingRules}, tree::AbstractGraph;
                          colorby::Symbol=:bonds,
                          maxpolys::Union{Integer,Nothing}=nothing,
                          title="Designable Search Tree",
                          filename=nothing,
                          colormap=:viridis,
                          figsize=nothing,
                          target=nothing)
    data = _target_tree_data(rules, tree, colorby)
    vmin, vmax = extrema(data.colorvalues)

    fig = Figure(; size=something(figsize, (900, _target_tree_height(data.nlevels))))
    _plot_target_tree!(fig.layout, data, tree; xmax=_target_tree_xmax(maxpolys, data.ps), vmin, vmax,
                       title, colormap, target)

    isnothing(filename) || save(filename, fig; backend=CairoMakie)
    return fig
end

"""
    plot_target_tree((rules, tree); kwargs...)

Plot the `(rules, tree)` pair returned by `create_designs_target` as is, see the method taking them apart.
"""
plot_target_tree((rules, tree)::Tuple{AbstractVector{<:BindingRules},AbstractGraph}; kwargs...) =
    plot_target_tree(rules, tree; kwargs...)

"""
    plot_target_trees(results, ncol; sharex=true, colorby=:bonds, maxpolys=nothing,
                      title="Designable Target Trees", subtitles=fill(nothing, length(results)),
                      filename=nothing, colormap=:viridis, figsize=nothing,
                      targets=fill(nothing, length(results)))

Plot each `(rules, tree)` pair in `results`, as returned by `create_designs_target`, as `plot_target_tree`
would, side by side in a grid of `ncol` columns filled row by row, under the overall `title`.

`subtitles` and `targets` hold the `title` and `target` of each plot, in the order of `results`; a
subtitle of `nothing` leaves its plot without one, and a `title` of `nothing` the figure.

With `sharex`, the axes of all plots run from 0 to `maxpolys`, by default the largest p of all trees, and
pan and zoom together, and only the plots with none below them label their x axis. Without it, each plot
has its own default `maxpolys`, as in `plot_target_tree`.

The colors mean the same in all plots, and a single colorbar to the right of the grid shows them.

Returns the figure, after saving it to `filename` with CairoMakie if given, as `plot_target_tree` does.
"""
function plot_target_trees(results::AbstractVector, ncol::Integer;
                           sharex::Bool=true,
                           colorby::Symbol=:bonds,
                           maxpolys::Union{Integer,Nothing}=nothing,
                           title="Designable Target Trees",
                           subtitles::AbstractVector=fill(nothing, length(results)),
                           filename=nothing,
                           colormap=:viridis,
                           figsize=nothing,
                           targets::AbstractVector=fill(nothing, length(results)))
    ncol >= 1 || throw(ArgumentError("ncol must be positive, got $ncol"))
    n = length(results)
    n >= 1 || throw(ArgumentError("results is empty"))
    length(subtitles) == n ||
        throw(ArgumentError("there are $(length(subtitles)) subtitles but $n (rules, tree) pairs"))
    length(targets) == n ||
        throw(ArgumentError("there are $(length(targets)) targets but $n (rules, tree) pairs"))

    data = [_target_tree_data(rules, tree, colorby) for (rules, tree) in results]
    vmin = minimum(d -> minimum(d.colorvalues), data)
    vmax = maximum(d -> maximum(d.colorvalues), data)
    allxmax = _target_tree_xmax(maxpolys, reduce(vcat, (d.ps for d in data)))

    # The rows of plots share the height in proportion to the height a single plot with the most levels in
    # each would get.
    nrow = cld(n, ncol)
    rowheights = [_target_tree_height(maximum(d -> d.nlevels, data[(r-1)*ncol+1:min(r*ncol, n)])) for r in 1:nrow]
    top = isnothing(title) ? 0 : 1
    fig = Figure(; size=something(figsize, (800ncol + 100, 40top + sum(rowheights))))
    isnothing(title) || Label(fig[1, 1:ncol], title; fontsize=24, font=:bold)

    axes = map(1:n) do k
        row, col = fldmod1(k, ncol)
        _plot_target_tree!(GridLayout(fig[top+row, col]), data[k], results[k][2];
                           xmax=sharex ? allxmax : _target_tree_xmax(maxpolys, data[k].ps), vmin, vmax,
                           title=subtitles[k], colormap, target=targets[k],
                           showx=!sharex || k + ncol > n, colorbar=false)
    end
    sharex && n > 1 && linkxaxes!(axes...)
    for (r, h) in enumerate(rowheights)
        rowsize!(fig.layout, top + r, Auto(false, h))
    end
    _integer_colorbar(fig[top+1:top+nrow, ncol+1], colormap, vmin, vmax, first(data).colorlabel)

    isnothing(filename) || save(filename, fig; backend=CairoMakie)
    return fig
end
