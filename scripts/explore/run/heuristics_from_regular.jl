# For tetromino pairs already searched with the full designability check by run/tetro_pairs.jl, the
# plain searches steered by `single` and by `exact` instead, each compared with the full one, and the
# output of the heuristics search of dtt_heuristics.jl rebuilt from all three.
#
# The rebuild replays the full search, so it computes only the checks the three searches leave open:
# the full check passed the rule sets the full search kept, and failed those it would have kept had
# they passed, and `single` passed every rule set the `single` and `exact` searches kept. It gives what
# `create_designs_target_heuristics` would, see `_search_designs_target_heuristics`, and is saved with
# those results in data/explore/heuristics-tetro-pairs.jld2. The two plain searches are saved in
# data/explore/regular-tetro-pairs.jld2, as `results[pair][:single]` and `results[pair][:exact]`.
# Searches already in either file are skipped.
#
# For every pair with a full search, without rebuilding anything: whether `single` and `exact` pass each
# rule set the full search kept, see `heuristic_finds`, which never runs the full check. Saved in
# data/explore/finds-tetro-pairs.jld2 as `results[pair] = (; single, exact)`, one `Bool` per rule set of
# the full search's `systems`, and printed per number of species. Pairs already there are not checked
# again.
#
# Run from the repository root with `julia --project=scripts -t 14 scripts/explore/run/heuristics_from_regular.jl`.

using Roly, Crafts, Clarabel, Convex, JLD2
include("../../lib/paths.jl")
include("../../lib/targetmatrices.jl")
include("../../lib/polyformfrommatrix.jl")
include("../dtt_heuristics.jl")

using Logging # Prevents warnings from optimizer
global_logger(ConsoleLogger(stderr, Logging.Error))

# As in run/tetro_pairs.jl, which made the full searches, and explore/run/heuristic_comparison.jl.
const PAIRS = (; maxsize=30, maxpolys=3000, greed=3)

sqtarget(group, name) = polyform_from_matrix(TargetInstantiators.SQMATS[group][name], UnitSquare)
loadresults(path, name) = isfile(path) ? load(path, name) : Dict{Any,Any}()

tetropairs = [("I", "J")]

full_results = load(datapath("tetro-pairs.jld2"), "tetro_results")
regular_path = datapath("explore", "regular-tetro-pairs.jld2")
regular = loadresults(regular_path, "results")
heuristics_path = datapath("explore", "heuristics-tetro-pairs.jld2")
heuristics = loadresults(heuristics_path, "results")

for pair in tetropairs
    targets = [sqtarget("tetrominoes", pair[1]), sqtarget("tetrominoes", pair[2])]
    full = full_results[pair]

    runs = get!(regular, pair, Dict{Symbol,Any}())
    for heuristic in (:single, :exact)
        haskey(runs, heuristic) && continue
        println(pair, ": ", heuristic)
        flush(stdout)
        runs[heuristic] = @time create_designs_target(targets, UnitSquare; PAIRS...,
                                                      single=heuristic == :single, exact=heuristic == :exact)
        jldsave(regular_path; results=regular, params=PAIRS)
    end

    for heuristic in (:single, :exact)
        same = sametree(runs[heuristic], full)
        println(pair, ": the ", heuristic, " search ", same ? "makes the same tree as" : "differs from", " the full one")
        same && continue
        println("  rule sets per number of species, full (a) against ", heuristic, " (b):")
        foreach(level -> println("  ", level), compare_trees(full[1], runs[heuristic][1]))
    end
    flush(stdout)

    haskey(heuristics, pair) && continue
    println(pair, ": rebuilding the heuristics search")
    flush(stdout)
    single_passes = Set(first(_canonicalforms(rules)) for heuristic in (:single, :exact)
                        for rules in runs[heuristic][1][2:end])
    counts = (exact=Threads.Atomic{Int}(0), single=Threads.Atomic{Int}(0), full=Threads.Atomic{Int}(0))
    heuristics[pair] = @time create_designs_target_heuristics(targets, UnitSquare; PAIRS..., known=full,
                                                              single_passes, counts)
    println(pair, ": ", length(heuristics[pair][3]), " rule sets checked, computing ",
            counts.exact[], " exact, ", counts.single[], " single and ", counts.full[], " full checks")
    flush(stdout)
    jldsave(heuristics_path; results=heuristics, params=PAIRS)
end

finds_path = datapath("explore", "finds-tetro-pairs.jld2")
finds = loadresults(finds_path, "results")

for pair in sort!(collect(keys(full_results)))
    if !haskey(finds, pair)
        targets = [sqtarget("tetrominoes", pair[1]), sqtarget("tetrominoes", pair[2])]
        finds[pair] = (; single=heuristic_finds(targets, full_results[pair]; maxsize=PAIRS.maxsize),
                       exact=heuristic_finds(targets, full_results[pair]; maxsize=PAIRS.maxsize, exact=true))
        jldsave(finds_path; results=finds, params=PAIRS)
    end

    systems = full_results[pair][1]
    (; single, exact) = finds[pair]
    println(pair, ": of ", length(systems), " rule sets kept by the full search, single passes ", count(single),
            " and exact ", count(exact))
    for m in sort!(unique(nspecies.(systems)))
        level = findall(rules -> nspecies(rules) == m, systems)
        println("  ", m, " species: ", length(level), " rule sets, single ", count(single[level]),
                ", exact ", count(exact[level]))
    end
    flush(stdout)
end
