using JuMP, CPLEX

file_path = "ExplForet_ampl2.dat"
if isfile(file_path)
    myFile = open(file_path)
    data = readlines(myFile)
    n = length(data)
    m = n
    T = zeros(n,m)
    for (i, line) in enumerate(data)
        line = split(line, " ")
        if length(line) > 0
            for (j, val) in enumerate(line)
                T[i,j] = parse(Float64, val)
            end
        end
    end
    close(myFile)
end

w1 = 2
w2 = 1
l = 3
g = 1.26157
m_quad = Model(CPLEX.Optimizer)
@variable(m_quad, x[1:n*m], Bin)

@objective(m_quad, Max, w1*sum(T[i, j] * (1 - x[(i - 1) * m + j]) for i in 1:n, j in 1:m) + w2 * g * l * sum(x[(i - 1) * m + j] * (sum(-x[(k - 1) * m + l] for k in 1:n, l in 1:m if (abs(k - i) + abs(l - j) == 1)) + 4) for i in 1:n, j in 1:m))

optimize!(m_quad)

feasibleSolutionFound = primal_status(m_quad) == MOI.FEASIBLE_POINT
isOptimal = termination_status(m_quad) == MOI.OPTIMAL
if feasibleSolutionFound
    vX = JuMP.value.(x)
    println("Valeur de l'objectif quad : ", JuMP.objective_value(m_quad))
    println("Node count: ", JuMP.node_count(m_quad))
else
    println("No feasible solution found")
end

m_lin = Model(CPLEX.Optimizer)
@variable(m_lin, 0 <= x[1:n*m] <= 1)
@variable(m_lin, y[1:n^2*m^2] >= 0)

for i in 1:n
    for j in 1:m
        for k in 1:n
            for l in 1:m
                if (abs(k - i) + abs(l - j) == 1)
                    @constraint(m_lin, x[(i - 1) * m + j] + x[(k - 1) * m + l] - y[(i - 1)*n*m^2 + (j - 1)n*m + (k - 1) * m + l] <= 1)
                else
                    @constraint(m_lin, y[(i - 1)*n*m^2 + (j - 1)n*m + (k - 1) * m + l] == 0)
                end
            end
        end
    end
end

@objective(m_lin, Max, w1*sum(T[i, j] * (1 - x[(i - 1) * m + j]) for i in 1:n, j in 1:m) + w2 * g * l * sum(4 * x[(i - 1) * m + j] - sum(y[(i - 1)*n*m^2 + (j - 1)n*m + (k - 1) * m + l] for k in 1:n, l in 1:m if abs(k - i) + abs(l - j) == 1) for i in 1:n, j in 1:m))

optimize!(m_lin)

feasibleSolutionFound = primal_status(m_lin) == MOI.FEASIBLE_POINT
isOptimal = termination_status(m_lin) == MOI.OPTIMAL
if feasibleSolutionFound
    vX = JuMP.value.(x)
    vY = JuMP.value.(y)
    println("Valeur de l'objectif : ", JuMP.objective_value(m_lin))
    println("Node count: ", JuMP.node_count(m_lin))
else
    println("No feasible solution found")
end


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