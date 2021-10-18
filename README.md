# Assignment

Routines to solve the linear assignment problem, and the k-best assignment problem.

This is a Julia port of the Matlab code for solving the 
[assignment problem](https://en.wikipedia.org/wiki/Assignment_problem) 
from the [TrackerComponentLibrary](https://github.com/USNavalResearchLaboratory/TrackerComponentLibrary) 
produced by D. F. Crouse [1], released into the public domain by the US Naval Research Laboratory.

This work is not affliated or endorsed by the Naval Research Laboratory.

### Examples
```julia
julia> using Assignment

julia> M=rand(1:100,3,4)
3×4 Matrix{Int64}:
 55  17  88  52
 29  38  90  20
 63  23  19  87

julia> sol = find_best_assignment(M)
AssignmentSolution(CartesianIndex.(1:3, [2, 4, 3]), 56)

julia> sum(M[sol])
56

julia> max_sol = find_best_assignment(M, true)
AssignmentSolution(CartesianIndex.(1:3, [1, 3, 4]), 232)

julia> sols = find_kbest_assigments(M, 5)
5-element Vector{Assignment.AssignmentSolution{Int64}}:
 AssignmentSolution(CartesianIndex.(1:4, [2, 4, 3, 1]), 56)
 AssignmentSolution(CartesianIndex.(1:4, [2, 1, 3, 4]), 65)
 AssignmentSolution(CartesianIndex.(1:4, [1, 4, 3, 2]), 94)
 AssignmentSolution(CartesianIndex.(1:4, [1, 4, 2, 3]), 98)
 AssignmentSolution(CartesianIndex.(1:4, [2, 4, 1, 3]), 100)
 
julia> max_sols = find_kbest_assigments(M, 5, true)
5-element Vector{Assignment.AssignmentSolution{Int64}}:
 AssignmentSolution(CartesianIndex.(1:4, [1, 3, 4, 2]), 232)
 AssignmentSolution(CartesianIndex.(1:4, [3, 2, 4, 1]), 213)
 AssignmentSolution(CartesianIndex.(1:4, [4, 3, 1, 2]), 205)
 AssignmentSolution(CartesianIndex.(1:4, [3, 1, 4, 2]), 204)
 AssignmentSolution(CartesianIndex.(1:4, [2, 3, 4, 1]), 194)
```

# References
[1] D. F. Crouse, "The Tracker Component Library: Free Routines for Rapid 
   Prototyping," IEEE Aerospace and Electronic Systems Magazine, vol. 32, 
   no. 5, pp. 18-27, May. 2017