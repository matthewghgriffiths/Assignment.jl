using Test
using Combinatorics
using Assignment

find_best(M, k) = size(M, 1) > size(M, 2) ? find_best(M', k) : partialsort([
    (sum(M[CartesianIndex.(axes(M, 1), p)]), p) for p in permutations(axes(M, 2), size(M, 1))], 1:k)

function test_kbest(M, k, sols)
    for ((cost, perm), sol) in zip(find_best(M, k), sols)
        @test cost ≈ sol.cost 
        @test perm == (sol.use_row ? sol.row4col : sol.col4row) 
    end
end

@testset "find best assignment" begin
      @testset "Float64" begin
          A = [ 0.891171  0.0320582   0.564188  0.8999    0.620615;
                0.166402  0.861136    0.201398  0.911772  0.0796335;
                0.77272   0.782759    0.905982  0.800239  0.297333;
                0.561423  0.170607    0.615941  0.960503  0.981906;
                0.748248  0.00799335  0.554215  0.745299  0.42637]

          sol = @inferred find_best_assignment(A)
          @test sol.row4col == [2, 3, 5, 1, 4]
          @test sol.cost ≈ 1.8375112
          @test sum(A[sol]) ≈ 1.8375112
      end

      @testset "Int64" begin 
          B = [ 24     1     8;
                 5     7    14;
                 6    13    20;
                12    19    21;
                18    25     2]

          sol = @inferred find_best_assignment(B)
          @test sol.row4col == [2, 1, 0, 0, 3]
          @test sol.col4row == [2, 1, 5]
          @test sol.cost == 8
          @test sum(B[sol]) == 8

          sol = @inferred find_best_assignment(B')
          @test sol.col4row == [2, 1, 0, 0, 3]
          @test sol.row4col == [2, 1, 5]
          @test sol.cost == 8
          @test sum(B'[sol]) == 8

          sol = @inferred find_best_assignment([ifelse(i-j==0,0,i*j) for i in 1:5, j in 1:5])
          @test sol.col4row == [1, 2, 3, 4, 5]
          @test sol.row4col == [1, 2, 3, 4, 5]
          @test sol.cost == 0
      end

      @testset "UInt16" begin
        M=UInt16[28092 44837 19882 39481 59139;
                 26039 46258 38932 51057 9;
                 11527 59487 61993 29072 8734;
                 10691 16977 12796 16370 14266;
                 5199  42319 34194 41332 16472]
        sol = @inferred find_best_assignment(M)
        @test sol.row4col == [3, 5, 4, 2, 1]
        @test sum(M[sol]) == 71139
        @test sol.cost    == 71139 % 2^16
    end
    
    @testset "UInt8" begin
        M=UInt8[67  228 135 197 244;
                112 44  84  206 31;
                225 103 231 225 227;
                170 37  135 9   130;
                110 22  133 77  96]

        sol = @inferred find_best_assignment(M)
        @test sol.row4col == [1, 5, 2, 4, 3]
        @test sum(M[sol]) == 343
        @test sol.cost   == 343 % 2^8
    end
end


@testset "find kbest assignments" begin
    @testset "Float64" begin
        A = [ 0.891171  0.0320582   0.564188  0.8999    0.620615;
              0.166402  0.861136    0.201398  0.911772  0.0796335;
              0.77272   0.782759    0.905982  0.800239  0.297333;
              0.561423  0.170607    0.615941  0.960503  0.981906;
              0.748248  0.00799335  0.554215  0.745299  0.42637]

        sols = @inferred find_kbest_assignments(A, 10)
        sol = sols[1]
        @test sol.row4col == [2, 3, 5, 1, 4]
        @test sol.cost ≈ 1.8375112
        @test sum(A[sol]) ≈ 1.8375112
        test_kbest(A, 10, sols)
    end

    @testset "Int64" begin 
        B = [ 24     1     8;
               5     7    14;
               6    13    20;
              12    19    21;
              18    25     2]

        sols = @inferred find_kbest_assignments(B, 10)
        sol = sols[1]
        @test sol.row4col == [2, 1, 0, 0, 3]
        @test sol.col4row == [2, 1, 5]
        @test sol.cost == 8
        @test sum(B[sol]) == 8
        test_kbest(B, 10, sort(sols, by=x -> (x.cost, x.col4row)))

        sols = @inferred find_kbest_assignments(B', 10)
        sol = sols[1]
        @test sol.col4row == [2, 1, 0, 0, 3]
        @test sol.row4col == [2, 1, 5]
        @test sol.cost == 8
        @test sum(B'[sol]) == 8
        test_kbest(B', 10, sort(sols, by=x -> (x.cost, x.row4col)))

        C = [ifelse(i-j==0,0,i*j) for i in 1:5, j in 1:5]
        sols = @inferred find_kbest_assignments(C, 10)
        sol = sols[1]
        @test sol.col4row == [1, 2, 3, 4, 5]
        @test sol.row4col == [1, 2, 3, 4, 5]
        @test sol.cost == 0
        test_kbest(C, 10, sort(sols, by=x -> (x.cost, x.row4col)))
    end

    @testset "UInt16" begin
      M=UInt16[28092 44837 19882 39481 59139;
               26039 46258 38932 51057 9;
               11527 59487 61993 29072 8734;
               10691 16977 12796 16370 14266;
               5199  42319 34194 41332 16472]
      sols = @inferred find_kbest_assignments(M, 10)
      sol = sols[1]
      @test sol.row4col == [3, 5, 4, 2, 1]
      @test sum(M[sol]) == 71139
      @test sol.cost    == 71139 % 2^16
  end
  
  @testset "UInt8" begin
      M=UInt8[67  228 135 197 244;
              112 44  84  206 31;
              225 103 231 225 227;
              170 37  135 9   130;
              110 22  133 77  96]

      sols = @inferred find_kbest_assignments(M, 10)
      sol = sols[1]
      @test sol.row4col == [1, 5, 2, 4, 3]
      @test sum(M[sol]) == 343
      @test sol.cost   == 343 % 2^8
  end
end



