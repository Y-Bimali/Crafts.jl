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

const PENTOMINOES = Dict( # Conway naming; "P1" is oriented like the letter P, "P2" is its mirror image
    "O" => [1 1 1 1 1],
    "P1" => [1 1; 1 1; 1 0], "P2" => [1 1; 1 1; 0 1],
    "Q1" => [1 0; 1 0; 1 0; 1 1], "Q2" => [0 1; 0 1; 0 1; 1 1], # Q1 is Golomb's L
    "R1" => [0 1 1; 1 1 0; 0 1 0], "R2" => [1 1 0; 0 1 1; 0 1 0], # R1 is the Game of Life R-pentomino
    "S1" => [0 0 1 1; 1 1 1 0], "S2" => [1 1 0 0; 0 1 1 1],
    "T" => [1 1 1; 0 1 0; 0 1 0],
    "U" => [1 0 1; 1 1 1],
    "V" => [1 0 0; 1 0 0; 1 1 1],
    "W" => [1 0 0; 1 1 0; 0 1 1],
    "X" => [0 1 0; 1 1 1; 0 1 0],
    "Y1" => [0 1; 1 1; 0 1; 0 1], "Y2" => [1 0; 1 1; 1 0; 1 0], # Y1 has its arm on the left
    "Z1" => [1 1 0; 0 1 0; 0 1 1], "Z2" => [0 1 1; 0 1 0; 1 1 0]
)

const SQMATS = Dict(
    "digits" => DIGITS,
    "tetrominoes" => TETROMINOES,
    "pentominoes" => PENTOMINOES
)

# ---------- Triangles
const TRMATS = Dict(
    "hexagon" => [1 1 1; 1 1 1]
)

end
