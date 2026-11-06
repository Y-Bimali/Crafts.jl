# Data for judging the `single` and `exact` heuristics: searches steered by the full designability
# check, see dtt_heuristics.jl, for the single tetrominoes, the single pentominoes, three tetromino
# pairs, and the eight shape. Each group is saved to its own file in data/explore after every target,
# and targets already in that file are skipped, so an interrupted run picks up where it stopped.
#
# Run from the repository root with `julia --project=scripts -t 13 scripts/explore/run/heuristic_comparison.jl`.

using Roly, Crafts, Clarabel, Convex, JLD2
include("../../lib/paths.jl")
include("../../lib/targetmatrices.jl")
include("../../lib/polyformfrommatrix.jl")
include("../dtt_heuristics.jl")

using Logging # Prevents warnings from optimizer
global_logger(ConsoleLogger(stderr, Logging.Error))

# Search parameters, as the results of create_designs_target these compare with were made: the singles
# and the eight shape without greed, the pairs as in run/tetro_pairs.jl.
const SINGLES = (; maxsize=30, maxpolys=3000, greed=nothing)
const PAIRS = (; maxsize=30, maxpolys=3000, greed=3)
const EIGHTSHAPE = (; maxsize=100, maxpolys=3000, greed=nothing)

sqtarget(group, name) = polyform_from_matrix(TargetInstantiators.SQMATS[group][name], UnitSquare)

"""
    run_heuristics!(file, jobs, params)

For each `key => targets` of `jobs` not yet in `file`, a file in data/explore, store
`create_designs_target_heuristics(targets, UnitSquare; params...)`, the triple
`(systems, network, notes)`, under `key`, saving the file after each. Returns the results so far.

The file holds `results`, a `Dict` of these triples, and `params`.
"""
function run_heuristics!(file, jobs, params)
    path = datapath("explore", file)
    results = isfile(path) ? load(path, "results") : Dict{Any,Any}()
    for (key, targets) in jobs
        haskey(results, key) && continue
        println(file, ": ", key)
        flush(stdout)
        results[key] = @time create_designs_target_heuristics(targets, UnitSquare; params...)
        jldsave(path; results, params)
    end
    return results
end

tetrominoes = sort!(collect(keys(TargetInstantiators.SQMATS["tetrominoes"])))
pentominoes = sort!(collect(keys(TargetInstantiators.SQMATS["pentominoes"])))
tetropairs = [("I", "O"), ("I", "S"), ("I", "T"), ("J", "S")]

# Cheapest first, and the eight shape, the one without a timing yet, last.
run_heuristics!("heuristics-tetrominoes.jld2", [name => [sqtarget("tetrominoes", name)] for name in tetrominoes], SINGLES)
run_heuristics!("heuristics-pentominoes.jld2", [name => [sqtarget("pentominoes", name)] for name in pentominoes], SINGLES)
run_heuristics!("heuristics-tetro-pairs.jld2",
                [(a, b) => [sqtarget("tetrominoes", a), sqtarget("tetrominoes", b)] for (a, b) in tetropairs], PAIRS)
run_heuristics!("heuristics-eightshape.jld2", ["8" => [sqtarget("digits", "8")]], EIGHTSHAPE)
