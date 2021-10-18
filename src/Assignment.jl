module Assignment

export find_best_assignment, find_kbest_assignments

import Base: <, >, <=, >=, ==, !=
using DataStructures

include("solution.jl")
include("solve.jl")
include("solve_kbest.jl")


end
