"""
find_kbest_assignments(C, k, maximize=false) -> kbest_sols

Find the k lowest [or highest] cost 2D assignments for the two-dimensional 
assignment problem with a rectangular cost matrix C.

# Example
```julia
julia> M=rand(1:100,3,4)
3×4 Matrix{Int64}:
 77  51  42  67
 72  53  47   4
 24  50  77  96

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

This is an implementation of Murty's method, which is described in [1].
The algorithm relies on the existence of a 2D assignment algorithm. The 2D
assignment algorithm of [2] is used. Additionally, the dual variable
inheritance methods described in [3] is used to reduce the computational
complexity of the technique.

Murty's algorithm runs 2D assignment algorithms a number of times with an
increasing number of constraints. Instances of MurtyData are stored in an 
ordered list implemented using the BinaryMinHeap class.

This code is a port of Matlab code released in the public domain by the
US Naval Research Laboratory [4]. 

This work is not affliated or endorsed by the Naval Research Laboratory.

# References
[1] K. G. Murty, "An algorithm for ranking all the assignments in order of
   increasing cost;" Operations Research; vol. 16; no. 3; pp. 682-687
   May-Jun. 1968.
[2] D. F. Crouse, "On Implementing 2D Rectangular Assignment Algorithms,"
   IEEE Transactions on Aerospace & Electronic Systems; vol. 52; no. 4
   pp. 1679-1696; Aug. 2016.
[3] M. L. Miller, H. S. Stone, & J. Cox, Ingemar, "Optimizing Murty's
   ranked assignment method;" IEEE Transactions on Aerospace &
   Electronic Systems; vol. 33; no. 3; pp. 851-862; Jul. 1997.
[4] D. F. Crouse, "The Tracker Component Library: Free Routines for Rapid 
   Prototyping," IEEE Aerospace and Electronic Systems Magazine, vol. 32, 
   no. 5, pp. 18-27, May. 2017

October 2013 David F. Crouse; Naval Research Laboratory; Washington D.C.
(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.
"""
function find_kbest_assignments(C::AbstractMatrix{T}, k::Integer, maximize=false) where T
    numRow, numCol = size(C)
	
	if numRow > numCol 
		sols = find_kbest_assignments(C', k, maximize)
		for (i, sol) in enumerate(sols)
			sols[i] = sol'
		end
		return sols
	end
    
	#The cost matrix must have all non-negative elements for the assignment
	#algorithm to work. This forces all of the elements to be positive. The
	#delta is added back in when computing the gain in the end.
    if maximize
        CDelta = maximum(C)
        C = -C.+ CDelta
	else 
        CDelta = minimum(C)
        C =  C .- CDelta
    end
	
	#Step 1: Find the best assignment. It will become the root from which()
	#everything else splits.
    #MUST PAD IT TO MAKE IT SQUARE; because we are using suboptimal
    #dual variables.
	numPad = numCol-numRow
	numPad > 0 && (C = [C; zeros(T, numPad, numCol)])
	
    LCHyp = MurtyData(C, numRow)
	

	#Check for feasibility.
    if isnothing(LCHyp.gainFull) 
        return AssignmentSolution{T, _dual_eltype(T)}[]
    end
	
	k_best = [AssignmentSolution(LCHyp, CDelta, maximize)]
	sizehint!(k_best, k)
	
	hyp_list = BinaryMinHeap([LCHyp])
	
	ScannedRows = falses(numCol)
	ScannedCol = falses(numCol)
	pred = zeros(Int, numCol)
	shortestPathCost = fill(typemax(float(T)), numCol)
	
    for curSweep in 2:k
        #We have to successively split the LC hypothesis either k times |
        #until splitting is no longer possible & get rid of the
        #hypothesis from the list.
		data = pop!(hyp_list)
		
		#Basically; it is all columns that have not been assigned in
		#one of the subproblems.
		col2Scan = @view data.col4rowLCFull[data.activeRow:end]
		
		#The loop ends at numVarRow; because we are not going to allow
		#variable hypotheses to be applied to dummy rows.
		for curRow=data.activeRow:data.numVarRow
			#forbiddenColumns contains a 1 for every column that one is 
			#not allowed to assign to curRow.
			if curRow==data.activeRow 
				forbiddenColumns = data.forbiddenActiveCol
			else
				forbiddenColumns = falses(numCol)
				forbiddenColumns[data.col4rowLCFull[curRow]] = true
			end

			#Get rid of the current assignment.
			row4colInit = data.row4colLCFull[:]
			col4rowInit = data.col4rowLCFull[:]
			row4colInit[col4rowInit[curRow]] = 0
			col4rowInit[curRow] = 0 

			splitHyp = ConstrainedShortestPath!(
				data.A, data.numVarRow, curRow, 
				forbiddenColumns, col4rowInit, row4colInit, col2Scan[:], data.u[:], data.v[:],
				ScannedRows, ScannedCol, pred, shortestPathCost
			)

			#If a valid assignment was produced.
			!isnothing(splitHyp.gainFull) && push!(hyp_list, splitHyp)

			#Remove the current assigned column from the list of
			#columns that can be scanned.
			col2Scan = @view col2Scan[col2Scan .!= data.col4rowLCFull[curRow]]
		end
		       	
        if !isempty(hyp_list)
			#Save the ordered best solutions.
			push!(k_best, AssignmentSolution(first(hyp_list), CDelta, maximize))
        else            
            break #All possible hypotheses have been enumerated.
        end
    end
    
	return k_best
end


struct MurtyData{T, M<:AbstractMatrix{T}, S}
	A::M
	numVarRow::Int
	activeRow::Int
	forbiddenActiveCol::BitVector
	col4rowLCFull::Vector{Int}
	row4colLCFull::Vector{Int}
	gainFull::Union{T, Nothing}
	u::Vector{S}
	v::Vector{S}
end

function MurtyData(A::AbstractMatrix, numVarRow::Integer)
	numRow, numCol = size(A)
	
	sol = find_best_assignment(A')
	forbiddenActiveCol = falses(numCol)
	if !isempty(sol.col4row) 
		forbiddenActiveCol[sol.col4row[1]] = 1
	end

	return MurtyData(
		A, numVarRow, 1, forbiddenActiveCol, 
		sol.col4row, sol.row4col, sol.cost, sol.v, sol.u)
end

Base.isless(a::MurtyData, b::MurtyData) = Base.isless(a.gainFull, b.gainFull)

for op in (:<, :>, :(<=), :(>=), :(==), :(!=))
	@eval begin
		$op(a::MurtyData, B::MurtyData) = $op(a.gainFull, b.gainFull)
		$op(a::MurtyData, B) = $op(a.gainFull, b)
		$op(a, B::MurtyData) = $op(a, b.gainFull)
	end
end


function AssignmentSolution(data::MurtyData{T, M, S}, CDelta=T(0), maximize=false) where {T, M, S}
	numCol = size(data.A, 2)
	numRow = data.numVarRow
	
	col4row = data.col4rowLCFull[1:numRow]
	row4col = data.row4colLCFull[:]
	if numRow < numCol
		for (i, j) in enumerate(row4col)
			if j > numRow
				row4col[i] = 0
			end
		end
	end
	gain =  maximize ? CDelta * numRow  - data.gainFull : CDelta * numRow + data.gainFull
	
		
	return AssignmentSolution{T, S}(col4row, row4col, T(gain), data.u[:], data.v[1:numRow])
end

function ConstrainedShortestPath!(
		C::AbstractMatrix{T}, 
		numVarRow::Integer, 
		activeRow::Integer, 
		forbiddenActiveCols::AbstractVector{Bool}, 
		col4row::AbstractVector{Int}, 
		row4col::AbstractVector{Int}, 
		col2Scan::AbstractVector{Int}, 
		u::AbstractVector{S}, 
		v::AbstractVector{S},
		ScannedRows::AbstractVector{Bool} = falses(size(C,1)),
		ScannedCol::AbstractVector{Bool} = falses(size(C,2)),
    	pred::AbstractVector{Int} = zeros(Int, size(C,2)),
		shortestPathCost::AbstractVector = fill(Inf, size(C,2))
	) where {T, S}
	
	numRow = size(C,1)
	numCol = size(C,2)
    numCol2Scan=length(col2Scan)
	ScannedRows[:] .= false
	ScannedCol[:] .= false
	pred[:] .= 0
	shortestPathCost[:] .= Inf
	
	sink = 0
	delta = S(0)
    curRow = activeRow
	closestColScan=0
	
	while sink==0
		ScannedRows[curRow]=1;
        
        # Scan all of the columns that have not already been scanned.
        minVal=Inf;
        @inbounds for curColScan=1:numCol2Scan
            curCol=col2Scan[curColScan];
			
			# Do not scan forbidden columns.
            curRow==activeRow && forbiddenActiveCols[curCol] && continue
            
            reducedCost=delta+C[curRow,curCol]-u[curRow]-v[curCol]
            if reducedCost < shortestPathCost[curCol]
                pred[curCol] = curRow
                shortestPathCost[curCol] = reducedCost
            end
            
            # Find the minimum unassigned column that was
            # scanned.
            if(shortestPathCost[curCol]<minVal)
                minVal=shortestPathCost[curCol];
                closestColScan=curColScan;
            end
        end
                
		# If the minimum cost column is not finite, then the problem is
		# not feasible.
		!isfinite(minVal) && return MurtyData(
			C, numVarRow, activeRow, forbiddenActiveCols, 
			Int[], Int[], nothing, u, v)
        
        closestCol=col2Scan[closestColScan];
        
        # Add the column to the list of scanned columns and delete it from
        # the list of columns to scan.		
        ScannedCol[closestCol] = 1
        numCol2Scan -= 1
		@inbounds for i in closestColScan:numCol2Scan
			col2Scan[i] = col2Scan[i + 1]
		end
        
        delta = S(shortestPathCost[closestCol])
        
        # If we have reached an unassigned row.
        if row4col[closestCol]==0 
            sink=closestCol
        else
            curRow=row4col[closestCol]
        end
    end
    
	# Dual Update Step
    # Update the first row in the augmenting path.
    u[activeRow] = u[activeRow] + delta;	
    # Update the rest of the rows in the agumenting path.
	@inbounds for (i, col) in enumerate(ScannedRows)
		if col != 0 && i != activeRow
			u[i] += delta - S(shortestPathCost[col4row[i]])
		end
	end
    # Update the scanned columns in the augmenting path.
	@inbounds for (i, row) in enumerate(ScannedCol)
		if row != 0
			v[i] -= delta - S(shortestPathCost[i])
		end
	end
	
	# Augmentation        
    # We have to remove node k from those that must be assigned.
    j = sink
    while true
        i = row4col[j] = pred[j]
		j, col4row[i] = col4row[i], j
		i==activeRow && break
    end
    
    # Calculate the gain that should be returned. Because the duality gap is
    # zero, the gain can be most easily obtained from the dual cost
    # function.
    gain=T(0)
	@inbounds for (i, j) in enumerate(col4row)
        gain=gain+C[i, j];
    end
	
	forbiddenActiveCols[col4row[activeRow]] = true
	
	return MurtyData(
		C, numVarRow, activeRow, forbiddenActiveCols, 
		col4row, row4col, gain, u, v)
end
