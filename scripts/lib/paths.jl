### Where the scripts keep their results: the data folder beside the Crafts.jl repository, outside it so
### that results never enter its history. Exploratory evaluations go in its `explore` subfolder. Setting
### the environment variable CRAFTS_DATA puts the folder elsewhere.

const DATA_DIR = get(ENV, "CRAFTS_DATA", normpath(joinpath(@__DIR__, "..", "..", "..", "data")))

"""
    datapath(parts...)

Return the path of `joinpath(parts...)` within `DATA_DIR`, creating the folder it lies in.
"""
function datapath(parts...)
    path = joinpath(DATA_DIR, parts...)
    mkpath(dirname(path))
    return path
end
