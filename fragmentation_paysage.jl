using JuMP, CPLEX

n = 30
couts = rand(1:10, n, n)

#couts = zeros(n, n)
# if isfile("Min_fragmentation_ampl.dat")
#     myFile = open("Min_fragmentation_ampl.dat")
#     data = readlines(myFile)
#     for i in 1:n
#         line = split(data[i], " ")
#         for j in 1:n
#             couts[i, j] = parse(Float64, line[j]) 
#         end
#     end
#     close(myFile)
# end

A_min = 100
A_max = 200
B = 1500

function distance(i, j, k, l)
    if i == k && j == l
        return 10000000
    end
    return sqrt((i - k) * (i - k) + (j - l) * (j - l))
end

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


m = Model(CPLEX.Optimizer)
@variable(m, x[1:n*n], Bin)
@variable(m, y[1:n*n, 1:n*n], Bin)

@constraint(m, A_min <= sum(x[i] for i in 1:n*n) <= A_max)
@constraint(m, sum(couts[i, j] * x[(i - 1) * n + j] * 10 for i in 1:n, j in 1:n) <= B)
@constraint(m, [i in 1:n*n, j in 1:n*n], y[i, j] <= x[j])
@constraint(m, [i in 1:n*n], sum(y[i, j] for j in 1:n*n) == x[i])

function Dinkelbach(lambda_0, epsilon)
    v = 15
    lambda = lambda_0
    vX = zeros(n*n)
    counter = 0
    while v < - epsilon || v > epsilon
        @objective(m, Min, sum(d[i, j] * y[i, j] for i in 1:n*n, j in 1:n*n) - lambda * sum(x[i] for i in 1:n*n))
        optimize!(m)
        vX = JuMP.value.(x)
        vY = JuMP.value.(y)
        v = JuMP.objective_value(m)
        println("Valeur optimale de P_lambda : ", v)
        lambda = sum(d[i, j] * vY[i, j] for i in 1:n*n for j in 1:n*n) / sum(vX[i] for i in 1:n*n)
        counter += 1
    end
    println("Valeur optimale  de P : ", lambda)
    println("Nombre d'itérations de l'algorithme : ", counter)
    return vX, lambda
end

X, lambda = Dinkelbach(0, 0.01)
println("Nombre de parcelles sélectionnées : ", sum(X[i] for i in 1:n*n))
for i in 1:n
    for j in 1:n
        print(X[(i - 1) * n + j])
        print(" | ")
    end
    println("")
end