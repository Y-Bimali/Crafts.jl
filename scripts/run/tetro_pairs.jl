# Short scripts that generate designable search trees results

using Roly, Crafts, Clarabel, Convex, JLD2
include("../lib/paths.jl")
include("../lib/targetmatrices.jl")
include("../lib/polyformfrommatrix.jl")
include("../lib/designabletargettree.jl")

using Logging # Prevents warnings from optimizer
global_logger(ConsoleLogger(stderr, Logging.Error))

# Flips the chirality of a binding rule set for UnitSquares
function chirality_swap(rules::BindingRules)
    swap(site) = site == 2 ? 4 : site == 4 ? 2 : site
    bonds = Roly.todict(rules)["bonds"]
    swapped = [[s1, swap(k1), s2, swap(k2)] for (s1, k1, s2, k2) in bonds]
    return BindingRules(swapped, UnitSquare)
end

# All pairs of tetrominoes
# - Equivalent pairs (i.e. differ only in chirality) are
#   skipped, rule graphs reversed, and populated.

distinct_pairs = [("I", "J"), ("I", "O"), ("I", "S"), ("I", "T"),
                  ("J", "L"), ("J", "O"), ("J", "S"), ("J", "T"), ("J", "Z"),
                  ("O", "S"), ("O", "T"),
                  ("S", "T"), ("S", "Z")]

chirality_matching = Dict(("I", "J") => ("I", "L"),
                          ("I", "S") => ("I", "Z"),
                          ("J", "O") => ("L", "O"),
                          ("J", "S") => ("L", "Z"),
                          ("J", "T") => ("L", "T"),
                          ("J", "Z") => ("L", "S"),
                          ("O", "S") => ("O", "Z"),
                          ("S", "T") => ("T", "Z"))

all_tetros = Dict(key => polyform_from_matrix(value, UnitSquare) for (key, value) in TargetInstantiators.SQMATS["tetrominoes"])

# tetro_results = Dict()
@load datapath("tetro-pairs.jld2") tetro_results

for pair in distinct_pairs
    if !haskey(tetro_results, pair)
        targets = [all_tetros[pair[1]], all_tetros[pair[2]]]
        println(pair[1] * pair[2])
        tetro_results[pair] = @time create_designs_target(targets, UnitSquare;
                                                        maxsize = 30, maxpolys = 3000,
                                                        greed = 3, verbose = true)
        @save datapath("tetro-pairs.jld2") tetro_results
    end
end

for (key, value) in chirality_matching
    tetro_results[value] = ([chirality_swap(rules) for rules in tetro_results[key][1]],
                            tetro_results[key][2])
end

@save datapath("tetro-pairs.jld2") tetro_results

# All digits (not greedy).
# - Equivalent targets (i.e. 2 and 5, 6 and 9) are skipped.
