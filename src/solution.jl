
"""
Stores the solution of an assigment problem

Can be directly used to index the appropriate values from the cost matrix 

The solution value is stored in the gain field
"""
struct AssignmentSolution{T}
	row4col::Vector{Int}
	col4row::Vector{Int}
	gain::Union{T, Nothing}
	v::Vector{T}
	u::Vector{T}
	use_row::Bool
end


function AssignmentSolution(row4col::Vector{Int}, col4row::Vector{Int}, gain, u::Vector{T}, v::Vector{T}) where {T}
	n = length(row4col)
	m = length(col4row)
	@assert n == length(v) 
	@assert m == length(u)
	AssignmentSolution{T}(row4col, col4row, gain, u, v, n <=m)
end

Base.adjoint(sol::AssignmentSolution) = AssignmentSolution(
	sol.col4row, sol.row4col, sol.gain, sol.u, sol.v, length(sol.col4row) <= length(sol.row4col))

Base.getindex(sol::AssignmentSolution, i::Integer) = sol.use_row ? CartesianIndex(i, sol.row4col[i]) : CartesianIndex(sol.col4row[i], i)

Base.eltype(sol::AssignmentSolution) = CartesianIndex{2}
Base.length(sol::AssignmentSolution) = sol.use_row ? length(sol.row4col) : length(sol.col4row)

Base.iterate(sol::AssignmentSolution) = Base.iterate(sol, 0)

function Base.iterate(sol::AssignmentSolution, i)
	i += 1
	i > length(sol) ? nothing : (sol[i], i)
end

Base.getindex(A::AbstractMatrix, sol::AssignmentSolution) = A[collect(sol)]

Base.show(io::IO, sol::AssignmentSolution) = if sol.use_row
    print(io, "AssignmentSolution(CartesianIndex.(1:$(length(sol.row4col)), $(sol.row4col)), $(sol.gain))")
else
    print(io, "AssignmentSolution(CartesianIndex.($(sol.col4row), 1:$(length(sol.col4row))), $(sol.gain))")
end
