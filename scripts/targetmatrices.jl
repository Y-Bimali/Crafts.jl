# a nested-Dict "filesystem" of instantiating matrices

module TargetInstantiators

export SQMATS, TRMATS # dicts for square and triangle formations

# ---------- Squares
const DIGITS = Dict(
    "0" => [1 1 1; 1 0 1; 1 0 1; 1 0 1; 1 1 1], #12
    "1" => [0 1 0; 1 1 0; 0 1 0; 0 1 0; 1 1 1], #8
    "2" => [1 1 1; 0 0 1; 1 1 1; 1 0 0; 1 1 1], #11, same as 5
    "3" => [1 1 1; 0 0 1; 1 1 1; 0 0 1; 1 1 1], #11
    "4" => [1 0 1; 1 0 1; 1 1 1; 0 0 1; 0 0 1], #9
    "5" => [1 1 1; 1 0 0; 1 1 1; 0 0 1; 1 1 1], #11, same as 2
    "6" => [1 1 1; 1 0 0; 1 1 1; 1 0 1; 1 1 1], #12, same as 9
    "7" => [1 1 1; 0 0 1; 0 0 1; 0 0 1; 0 0 1], #7
    "8" => [1 1 1; 1 0 1; 1 1 1; 1 0 1; 1 1 1], #13
    "9" => [1 1 1; 1 0 1; 1 1 1; 0 0 1; 1 1 1], #12, same as 6
)

const TETROMINOES = Dict(
    "I" => [1 1 1 1],
    "O" => [1 1; 1 1],
    "T" => [1 1 1; 0 1 0],
    "J" => [1 1 1; 0 0 1],
    "L" => [0 0 1; 1 1 1],
    "S" => [0 1 1; 1 1 0],
    "Z" => [1 1 0; 0 1 1]
)

const SQMATS = Dict(
    "digits" => DIGITS,
    "tetrominoes" => TETROMINOES
)

# ---------- Triangles
const TRMATS = Dict(
    "hexagon" => [1 1 1; 1 1 1]
)

end
