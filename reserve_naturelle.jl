using JuMP, Gurobi#, AmplNLWriter
using Distributions

#Params
sz_map = 6
n_species = 12
unif_rate = 0.9

## Generating the terrain
is_generated = zeros(sz_map, sz_map)
is_seen = zeros(sz_map, sz_map)
next_to_price = Array([(rand(0:sz_map), rand(0:sz_map))])
price = rand(1:10, sz_map, sz_map) #couts

while length(next_to_price) > 0
	next_tile = popfirst!(next_to_price)
	i = next_tile[1]
	j = next_tile[2]
	is_generated[i,j] = 1
	n_generated_neighs = 0
	tot_neigh_price = 0
	## Make sure we check the neighbours later
	if i > 1 && is_seen[i-1, j] == 0
		is_seen[i-1, j] = 1
		push!(next_to_price, (i-1, j))
	elseif i > 1 && is_generated[i-1, j] == 1
		n_generated_neighs += 1
		tot_neigh_price += price[i-1,j]
	end
	if i < sz_map - 1 && is_seen[i+1, j] == 0
		is_seen[i+1, j] = 1
		push!(next_to_price, (i+1, j))
	elseif i < sz_map - 1 && is_generated[i+1, j] == 1
		n_generated_neighs += 1
		tot_neigh_price += price[i+1,j]
	end
	if j > 1 && is_seen[i, j-1] == 0
		is_seen[i, j-1] = 1
		push!(next_to_price, (i, j-1))
	elseif j > 1 && is_generated[i, j-1] == 1
		n_generated_neighs += 1
		tot_neigh_price += price[i,j-1]
	end
	if j < sz_map - 1 && is_seen[i, j+1] == 0
		is_seen[i, j+1] = 1
		push!(next_to_price, (i, j+1))
	elseif j < sz_map - 1 && is_generated[i, j+1] == 1
		n_generated_neighs += 1
		tot_neigh_price += price[i,j+1]
	end
	gamma = unif_rate * n_generated_neighs / 3
	price[i,j] = trunc(Int, (1 - gamma) * price[i,j] + gamma * tot_neigh_price / max(1, n_generated_neighs))
end

## Generating survival
species_good_habitat_rate = rand(Uniform(0.0,0.2), n_species) #alpha
zones = zeros(n_species, sz_map, sz_map)
habitability = rand(Uniform(0,1), sz_map, sz_map)
for spe in 1:n_species
	for i in 1:sz_map
        for j in 1:sz_map
            if habitability[i,j] < rand(Uniform(0,1))/20 + species_good_habitat_rate[spe]
                zones[spe,i,j] = 2*habitability[i,j] + rand(Uniform(0,1))/5
            end
        end
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

m = Model(Gurobi.Optimizer)
@variable(m, x[1:n*n], Bin)
@variable(m, q[1:n*n], Bin)

@constraint(m, [i in 1:3], log(1 - alpha[i]) >= sum(log(1 - proba[i, k, l]) * q[(k - 1) * n + l] for k in 1:n, l in 1:n))
@constraint(m, [i in 4:6], log(1 - alpha[i]) >= sum(log(1 - proba[i, k, l]) * x[(k - 1) * n + l] for k in 1:n, l in 1:n))
@constraint(m, [k in 2:n-1, l in 2:n-1], 9 * q[(k - 1) * n + l] <= sum(x[(i - 1) * n + j] for i in k-1:k+1, j in l-1:l+1))
@constraint(m, [l in 1:n], q[l] == 0)
@constraint(m, [l in 1:n], q[(n - 1) * n + l] == 0)
@constraint(m, [k in 1:n], q[(k - 1) * n + 1] == 0)
@constraint(m, [k in 1:n], q[(k - 1) * n + n] == 0)

@objective(m, Min, sum(couts[i][j] * x[(i - 1) * n + j] for i in 1:n, j in 1:n))

optimize!(m)

feasibleSolutionFound = primal_status(m) == MOI.FEASIBLE_POINT
isOptimal = termination_status(m) == MOI.OPTIMAL
if feasibleSolutionFound
    vX = JuMP.value.(x)
    vQ = JuMP.value.(q)
    println("Valeur de l'objectif : ", JuMP.objective_value(m))
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