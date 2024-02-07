using JuMP, CPLEX

function run(n, p_d, p_c, proba, alpha)
    m = Model(CPLEX.Optimizer)
    @variable(m, x[1:n*n], Bin)
    @variable(m, q[1:n*n], Bin)

    @constraint(m, [i in 1:p_d], log(1 - alpha[i]) >= sum(log(1 - proba[i, k, l]) * q[(k - 1) * n + l] for k in 1:n, l in 1:n))
    @constraint(m, [i in p_d+1:p_d+p_c], log(1 - alpha[i]) >= sum(log(1 - proba[i, k, l]) * x[(k - 1) * n + l] for k in 1:n, l in 1:n))
    @constraint(m, [k in 2:n-1, l in 2:n-1], 9 * q[(k - 1) * n + l] <= sum(x[(i - 1) * n + j] for i in k-1:k+1, j in l-1:l+1))
    @constraint(m, [l in 1:n], q[l] == 0)
    @constraint(m, [l in 1:n], q[(n - 1) * n + l] == 0)
    @constraint(m, [k in 1:n], q[(k - 1) * n + 1] == 0)
    @constraint(m, [k in 1:n], q[(k - 1) * n + n] == 0)

    @objective(m, Min, sum(couts[i][j] * x[(i - 1) * n + j] for i in 1:n, j in 1:n))

    start_time = time()
    optimize!(m)
    end_time = time()

    feasibleSolutionFound = primal_status(m) == MOI.FEASIBLE_POINT
    if feasibleSolutionFound
        vX = JuMP.value.(x)
        vQ = JuMP.value.(q)
        v = JuMP.objective_value(m)
        # println("Valeur de l'objectif : ", JuMP.objective_value(m))
        proba_survie = zeros(p_d + p_c)
        for i in 1:p_d
            proba_survie[i] = 1 - prod(1 - proba[i, k, l] * vQ[(k - 1) * n + l] for k in 1:n, l in 1:n)
        end
        for i in p_d+1:p_d+p_c
            proba_survie[i] = 1 - prod(1 - proba[i, k, l] * vX[(k - 1) * n + l] for k in 1:n, l in 1:n)
        end
        return end_time - start_time, MOI.get(m, MOI.NodeCount()), JuMP.objective_value(m), proba_survie
    else
        # println("No feasible solution found")
        return false
    end 
end

n = 10
p = 6
proba = zeros(p, n, n)

if isfile("RODD-probabilites.txt")
    myFile = open("RODD-probabilites.txt")
    data = readlines(myFile)
    for line in data
        line = split(line, " ")
        if length(line) > 0
            proba[parse(Int64, line[1]), parse(Int64, line[2]), parse(Int64, line[3])] = parse(Float64, line[4]) 
        end
    end
    close(myFile)
end

couts = [[6, 6, 6, 4, 4, 4, 4, 8, 8, 8],
         [6, 6, 6, 4, 4, 4, 4, 8, 8, 8],
         [6, 6, 6, 4, 4, 4, 4, 8, 8, 8],
         [5, 5, 5, 3, 3, 3, 3, 7, 7, 7],
         [5, 5, 5, 3, 3, 3, 3, 7, 7, 7],
         [5, 5, 5, 3, 3, 3, 3, 7, 7, 7],
         [5, 5, 5, 3, 3, 3, 3, 7, 7, 7],
         [4, 4, 4, 6, 6, 6, 6, 5, 5, 5],
         [4, 4, 4, 6, 6, 6, 6, 5, 5, 5],
         [4, 4, 4, 6, 6, 6, 6, 5, 5, 5]]

alpha = [0.8, 0.8, 0.8, 0.6, 0.6, 0.6]

run(n, 3, 3, proba, alpha)



# C = 300
# mlp = Model(() -> AmplNLWriter.Optimizer("bonmin"))
# @variable(mlp, xlp[1:n*n], Bin)
# @variable(mlp, qlp[1:n*n], Bin)

# @constraint(mlp, [k in 2:n-1, l in 2:n-1], 9 * qlp[(k - 1) * n + l] <= sum(xlp[(i - 1) * n + j] for i in k-1:k+1, j in l-1:l+1))
# @constraint(mlp, [l in 1:n], qlp[l] == 0)
# @constraint(mlp, [l in 1:n], qlp[(n - 1) * n + l] == 0)
# @constraint(mlp, [k in 1:n], qlp[(k - 1) * n + 1] == 0)
# @constraint(mlp, [k in 1:n], qlp[(k - 1) * n + n] == 0)
# @constraint(mlp, sum(couts[i][j] * xlp[(i - 1) * n + j] for i in 1:n, j in 1:n) <= C)

# function define_expr(x_var, q_var)
#     expr = 0
#     for i in 1:1
#         prod = 1
#         for k in 1:5
#             for l in 1:n
#                 if proba[i, k, l] != 0.0
#                     prod = prod * (1 - proba[i, k, l] * q_var[(k - 1) * n + l])
#                 end
#             end
#         end
#         expr += prod
#     end
#     return expr
# end

# println(define_expr(xlp, qlp))
# @objective(mlp, Max, define_expr(xlp, qlp))

# optimize!(mlp)

# feasibleSolutionFoundLp = primal_status(mlp) == MOI.FEASIBLE_POINT
# isOptimalLp = termination_status(mlp) == MOI.OPTIMAL
# if feasibleSolutionFoundLp
#     vXLp = JuMP.value.(xlp)
#     vQLp = JuMP.value.(qlp)
#     println("Valeur de l'objectif : ", JuMP.objective_value(mlp))
# else
#     println("No feasible solution found huh")
# end