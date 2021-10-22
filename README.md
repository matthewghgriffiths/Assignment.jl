# Assignment

Routines to solve the linear assignment problem, and the _k_-best assignment problem.

A variant of the Jonker-Volgenant algorithm is used to solve the linear assignment problem, this algorithm is in general faster than the Hungarian/Kuhn-Munkres algorithm. Murty's algorithm combined with the modified Jonker-Volgenant algorithm is used to solve the _k_-best assignment problem. Algorithms are described in detail in [1].

This is a Julia port of the Matlab code for solving the 
[assignment problem](https://en.wikipedia.org/wiki/Assignment_problem) 
from the [TrackerComponentLibrary](https://github.com/USNavalResearchLaboratory/TrackerComponentLibrary) 
produced by D. F. Crouse [2], released into the public domain by the US Naval Research Laboratory.

This work is not affliated with or endorsed by the US Naval Research Laboratory.

### Examples
```julia
julia> using Assignment

julia> M=rand(1:100,3,4)
3×4 Matrix{Int64}:
 77  51  42  67
 72  53  47   4
 24  50  77  96

julia> sol = find_best_assignment(M)
AssignmentSolution(CartesianIndex.(1:3, [3, 4, 1]), 70)

julia> sum(M[sol])
70

julia> max_sol = find_best_assignment(M', true)
AssignmentSolution(CartesianIndex.([1, 2, 4], 1:3), 226)

julia> sols = find_kbest_assigments(M, 5)
5-element Vector{Assignment.AssignmentSolution{Int64, Int64}}:
 AssignmentSolution(CartesianIndex.(1:3, [3, 4, 1]), 70)
 AssignmentSolution(CartesianIndex.(1:3, [2, 4, 1]), 79)
 AssignmentSolution(CartesianIndex.(1:3, [3, 4, 2]), 96)
 AssignmentSolution(CartesianIndex.(1:3, [3, 2, 1]), 119)
 AssignmentSolution(CartesianIndex.(1:3, [2, 3, 1]), 122)
 
julia> max_sols = find_kbest_assignments(M, 5, true)
5-element Vector{Assignment.AssignmentSolution{Int64, Int64}}:
 AssignmentSolution(CartesianIndex.(1:3, [1, 2, 4]), 226)
 AssignmentSolution(CartesianIndex.(1:3, [1, 3, 4]), 220)
 AssignmentSolution(CartesianIndex.(1:3, [2, 1, 4]), 219)
 AssignmentSolution(CartesianIndex.(1:3, [4, 1, 3]), 216)
 AssignmentSolution(CartesianIndex.(1:3, [3, 1, 4]), 210)
```

# References
[1] D. F. Crouse, "Advances in displaying uncertain estimates of multiple
   targets," in Proceedings of SPIE: Signal Processing, Sensor Fusion, and
   Target Recognition XXII, vol. 8745, Baltimore, MD, Apr. 2013.
[2] D. F. Crouse, "The Tracker Component Library: Free Routines for Rapid 
   Prototyping," IEEE Aerospace and Electronic Systems Magazine, vol. 32, 
   no. 5, pp. 18-27, May. 2017