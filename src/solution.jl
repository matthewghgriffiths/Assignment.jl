"""
Stores the solution of an assigment problem

Can be directly used to index the appropriate values from the cost matrix 
"""
struct AssignmentSolution{T, S}
	row4col::Vector{Int}
	col4row::Vector{Int}
	cost::Union{T, Nothing}
	v::Vector{S}
	u::Vector{S}
	use_row::Bool

    # AssignmentSolution(row4col, col4row, cost::Nothing, u::Vector{S}, v::Vector{S}, use_row) where S = new{S, S}(row4col, col4row, cost, u, v, use_row)
    AssignmentSolution{T, S}(row4col, col4row, cost::T, u::Vector{S}, v::Vector{S}, use_row) where {T, S} = new{T, S}(row4col, col4row, cost, u, v, use_row)
    AssignmentSolution{T, S}(row4col, col4row, cost::Nothing, u::Vector{S}, v::Vector{S}, use_row) where {T, S} = new{T, S}(row4col, col4row, cost, u, v, use_row)
    AssignmentSolution(row4col, col4row, cost::T, u::Vector{S}, v::Vector{S}, use_row) where {T, S} = new{T, S}(row4col, col4row, cost, u, v, use_row)
end

function AssignmentSolution{T, S}(row4col, col4row, cost::Union{T, Nothing}, u::Vector{S}, v::Vector{S}) where {T, S}
	n = length(row4col)
	m = length(col4row)
	@assert n == length(v) 
	@assert m == length(u)
	AssignmentSolution{T, S}(row4col, col4row, cost, u, v, n <=m)
end

AssignmentSolution(row4col::Vector{Int}, col4row::Vector{Int}, cost::T, u::Vector{S}, v::Vector{S}) where {T, S} = 
    AssignmentSolution{isnothing(cost) ? S : T, S}(row4col, col4row, cost, u, v)

Base.adjoint(sol::AssignmentSolution{T, S}) where {T, S} = AssignmentSolution{T, S}(
	sol.col4row, sol.row4col, sol.cost, sol.u, sol.v, length(sol.col4row) <= length(sol.row4col))

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
    print(io, "AssignmentSolution(CartesianIndex.(1:$(length(sol.u)), $(sol.row4col)), $(sol.cost))")
else
    print(io, "AssignmentSolution(CartesianIndex.($(sol.col4row), 1:$(length(sol.v))), $(sol.cost))")
end
