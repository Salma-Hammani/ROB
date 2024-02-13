using JuMP, CPLEX

function distance(i, j, k, l)
    if i == k && j == l
        return 10000000
    end
    return sqrt((i - k) * (i - k) + (j - l) * (j - l))
end

function run(n, couts, A_min, A_max, B)
    d = zeros(n*n, n*n)
    for i in 1:n
        for j in 1:n
            for k in 1:n
                for l in 1:n
                    d[(i - 1) * n + j, (k - 1) * n + l] = distance(i, j, k, l)
                end
            end
        end
    end

    m = Model(optimizer_with_attributes(CPLEX.Optimizer, "CPXPARAM_ScreenOutput" => 0))
    @variable(m, x[1:n*n], Bin)
    @variable(m, y[1:n*n, 1:n*n], Bin)

    @constraint(m, A_min <= sum(x[i] for i in 1:n*n) <= A_max)
    @constraint(m, sum(couts[i, j] * x[(i - 1) * n + j] * 10 for i in 1:n, j in 1:n) <= B)
    @constraint(m, [i in 1:n*n, j in 1:n*n], y[i, j] <= x[j])
    @constraint(m, [i in 1:n*n], sum(y[i, j] for j in 1:n*n) == x[i])

    v = 15
    lambda = 0.
    vX = zeros(n*n)
    counter = 0
    print(string(counter) * ".. ")
    t = 0
    nodes = 0
    while v < - 0.01 || v > 0.01
        @objective(m, Min, sum(d[i, j] * y[i, j] for i in 1:n*n, j in 1:n*n) - lambda * sum(x[i] for i in 1:n*n))
        start_time = time()
        optimize!(m)
        end_time = time()
        feasibleSolutionFound = primal_status(m) == MOI.FEASIBLE_POINT
        if feasibleSolutionFound
            vX = JuMP.value.(x)
            vY = JuMP.value.(y)
            v = JuMP.objective_value(m)
            # println("Valeur optimale de P_lambda : ", v)
            nodes += MOI.get(m, MOI.NodeCount())
            t += end_time - start_time
            lambda = sum(d[i, j] * vY[i, j] for i in 1:n*n for j in 1:n*n) / sum(vX[i] for i in 1:n*n)
            counter += 1
            print(string(counter) * ".. ")
        else
            println("")
            return false
        end
    end
    # println("Valeur optimale  de P : ", lambda)
    # println("Nombre d'itérations de l'algorithme : ", counter)
    println("")
    return t, nodes, counter, lambda, vX
end

# couts = zeros(n, n)
# if isfile("Min_fragmentation_ampl.dat")
#     myFile = open("Min_fragmentation_ampl.dat")
#     data = readlines(myFile)
#     for i in 1:n
#         line = split(data[i], " ")
#         for j in 1:n
#             couts[i, j] = parse(Float64, line[j + 1]) 
#         end
#     end
#     close(myFile)
# end
# time, nodes, counter, lambda, X = run(n, couts, A_min, A_max, B)
# println("Nombre de parcelles sélectionnées : ", sum(X[i] for i in 1:n*n))
# for i in 1:n
#     for j in 1:n
#         print(X[(i - 1) * n + j])
#         print(" | ")
#     end
#     println("")
# end

for n in 38:40
    println("n = " * string(n))
    couts = rand(1:10, n, n)
    A_min = rand(n * 2 : n * 7)
    A_max = rand(A_min+1:A_min + 5)
    B = rand(A_min : 30 * A_min)
    result = run(n , couts, A_min, A_max, B)
    while (result == false)
        couts = rand(1:10, n, n)
        A_min = rand(n * 2 : n * 7)
        A_max = rand(A_min+1:A_min + 5)
        B = rand(A_min : 30 * A_min)
        result = run(n , couts, A_min, A_max, B)
    end
    file = open("fragmentation.csv", "a")
    write(file, string(floor(n)))
    write(file, ";")
    write(file, string(floor(A_min)))
    write(file, ";")
    write(file, string(floor(A_max)))
    write(file, ";")
    write(file, string(floor(B)))
    write(file, ";")
    write(file, string(round(result[1], digits=3)))
    write(file, ";")
    write(file, string(floor(result[2])))
    write(file, ";")
    write(file, string(floor(result[3])))
    write(file, ";")
    write(file, string(round(result[4], digits=3)))
    write(file, "\n")
    close(file)
end

