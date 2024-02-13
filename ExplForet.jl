using JuMP, CPLEX

function runQuad(n, m, T, w1, w2, g, l, lower_limit=0)
    m_quad = Model(optimizer_with_attributes(CPLEX.Optimizer, "CPXPARAM_ScreenOutput" => 0))
    @variable(m_quad, x[1:n*m], Bin)

    @objective(m_quad, Max, w1*sum(T[i, j] * (1 - x[(i - 1) * m + j]) for i in 1:n, j in 1:m) + w2 * g * l * sum(x[(i - 1) * m + j] * (sum(-x[(k - 1) * m + l] for k in 1:n, l in 1:m if (abs(k - i) + abs(l - j) == 1)) + 4) for i in 1:n, j in 1:m))

    if lower_limit != 0
        @constraint(m_quad, sum(x[k] for k in 1:n*m) >= 60)
    end
    start_time = time()
    optimize!(m_quad)
    end_time = time()

    feasibleSolutionFound = primal_status(m_quad) == MOI.FEASIBLE_POINT
    if feasibleSolutionFound
        vX = JuMP.value.(x)
        v = JuMP.objective_value(m_quad)
        noeuds = JuMP.node_count(m_quad)
        return end_time - start_time, noeuds, v, vX
    else
        return false
    end
end

function runLin(n, m, T, w1, w2, g, l, lower_limit=0)
    m_lin = Model(optimizer_with_attributes(CPLEX.Optimizer, "CPXPARAM_ScreenOutput" => 0))
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
    
    if lower_limit != 0
        @constraint(m_lin, sum(x[k] for k in 1:n*m) >= 60)
    end

    @objective(m_lin, Max, w1*sum(T[i, j] * (1 - x[(i - 1) * m + j]) for i in 1:n, j in 1:m) + w2 * g * l * sum(4 * x[(i - 1) * m + j] - sum(y[(i - 1)*n*m^2 + (j - 1)n*m + (k - 1) * m + l] for k in 1:n, l in 1:m if abs(k - i) + abs(l - j) == 1) for i in 1:n, j in 1:m))
    
    start_time = time()
    optimize!(m_lin)
    end_time = time()

    feasibleSolutionFound = primal_status(m_lin) == MOI.FEASIBLE_POINT
    if feasibleSolutionFound
        vX = JuMP.value.(x)
        vY = JuMP.value.(y)
        # println("Valeur de l'objectif : ", JuMP.objective_value(m_lin))
        # println("Node count: ", JuMP.node_count(m_lin))
        v = JuMP.objective_value(m_lin)
        noeuds = JuMP.node_count(m_lin)
        return end_time - start_time, noeuds, v, vX, vY
    else
        return false
    end
end


# file_path = "ExplForet_ampl2.dat"
# if isfile(file_path)
#     myFile = open(file_path)
#     data = readlines(myFile)
#     n = length(data)
#     m = n
#     T = zeros(n,m)
#     for (i, line) in enumerate(data)
#         line = split(line, " ")
#         if length(line) > 0
#             for (j, val) in enumerate(line)
#                 T[i,j] = parse(Float64, val)
#             end
#         end
#     end
#     close(myFile)
# end

# w1 = 5
# w2 = 2
# l = 3
# g = 1.26157
# lb = 60

# quad_time, quad_noeuds, quad_val, quad_x = runQuad(n, m, T, w1, w2, g, lb)
# println(quad_val)
# lin_time, lin_noeuds, lin_val, lin_x, lin_y = runLin(n, m, T, w1, w2, g, lb)
# println(lin_val)

# for i in 1:n
#     for j in 1:n
#         print(quad_x[(i - 1) * n + j])
#         print(" | ")
#     end
#     println("")
# end

for n in 35:50
    println("n = " * string(n))
    m = n
    T = rand(50:100, n, m)
    w1 = rand(1:5)
    w2 = rand(1:5)
    l = 3
    g = 1.26157
    lb = rand(floor(n*m/4):3*floor(n*m/4))
    result1 = runQuad(n, m, T, w1, w2, g, l, lb)
    result2 = runLin(n, m, T, w1, w2, g, l, lb)
    while (result1 == false)
        T = rand(50:100, n, m)
        w1 = rand(1:5)
        w2 = rand(1:5)
        l = 3
        g = 1.26157
        lb = rand(floor(n*m/4):3*floor(n*m/4))
        result1 = runQuad(n, m, T, w1, w2, g, l, lb)
        result2 = runLin(n, m, T, w1, w2, g, l, lb)
    end
    quad_time, quad_noeuds, quad_val, quad_x = result1
    lin_time, lin_noeuds, lin_val, lin_x, lin_y = result2

    file = open("epxloitationForet.csv", "a")
    write(file, string(floor(n)))
    write(file, ";")
    write(file, string(floor(m)))
    write(file, ";")
    write(file, string(floor(w1)))
    write(file, ";")
    write(file, string(floor(w2)))
    write(file, ";")
    write(file, string(floor(lb)))
    write(file, ";")
    write(file, string(round(g, digits=2)))
    write(file, ";")
    write(file, string(round(quad_time, digits=2)))
    write(file, ";")
    write(file, string(floor(quad_noeuds)))
    write(file, ";")
    write(file, string(round(quad_val, digits=2)))
    write(file, ";")
    write(file, string(round(lin_time, digits=2)))
    write(file, ";")
    write(file, string(floor(lin_noeuds)))
    write(file, ";")
    write(file, string(round(lin_val, digits=2)))
    write(file, "\n")
    close(file)
end
