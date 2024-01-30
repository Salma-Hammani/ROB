using JuMP, CPLEX

P = 5000 # nombre d'individu de même sexe
M = 20 # nombre de gènes par chromosome
K = 2 # nombre de chromosomes par individu
N = 15 # nombre d'allèles par gène

nMat = zeros(2 * P, M, N)
for k in 1:2*P
    for i in 1:M
        nMat[k, i, 1] = rand(0:K)
        for j in 2:(N-1)
            nMat[k, i, j] = rand(0:K - sum(nMat[k, i, l] for l in 1:j-1))
        end
        nMat[k, i, N] = K - sum(nMat[k, i, l] for l in 1:N-1)
    end
end

# if isfile("DivGenetique_ampl.dat")
#     myFile = open("DivGenetique_ampl.dat")
#     data = readlines(myFile)
#     for i in 1:2*M*2*P
#         line = split(data[i], " ")
#         individu = parse(Int64, line[1])
#         gene = parse(Int64, line[3])
#         allele = parse(Int64, line[5])
#         nMat[individu, gene, allele] += 1
#     end
#     close(myFile)
# end

h = 100
theta1 = 0.001
theta = zeros(h)
for r in 1:h
    theta[r] = theta1 ^ ((h - r) / (h - 1))
end

m = Model(CPLEX.Optimizer)
@variable(m, d[1:2*P], Int)
@variable(m, 0 <= t[1:M, 1:N] <= 1)
@variable(m, 0 <= z[1:M, 1:N] <= 1)

@constraint(m, sum(d[k] for k in 1:P) == sum(d[k] for k in P+1:2*P))
@constraint(m, sum(d[k] for k in 1:P) == 2 * P)
@constraint(m, [r in 1:h, i in 1:M, j in 1:N], log(theta[r]) + (z[i, j] - theta[r]) / theta[r] >= sum(d[k] * log(1/2) for k in 1:2*P if nMat[k, i, j] == 1))
@constraint(m, [i in 1:M, j in 1:N], t[i, j] >= z[i, j] - sum(d[k] for k in 1:2*P if nMat[k, i, j] == 2))
#@constraint(m, [k in 1:2*P], d[k] <= 3)
@constraint(m, [k in 1:2*P], d[k] >= 0)

@objective(m, Min, sum(t[i, j] for i in 1:M, j in 1:N))

optimize!(m)
feasibleSolutionFound = primal_status(m) == MOI.FEASIBLE_POINT
isOptimal = termination_status(m) == MOI.OPTIMAL
if feasibleSolutionFound
    vD = JuMP.value.(d)
    println("Valeur de la borne inf: ", JuMP.objective_value(m))
    println("Nombre de noeuds explorés ", MOI.get(m, MOI.NodeCount()))
else
    println("No feasible solution found")
end

