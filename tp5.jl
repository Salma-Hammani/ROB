using JuMP, Gurobi

#Parameters
T = 10
SURF = 1
P = 40
a_bar = 2
l_bar = 2

# Instance
crops = Dict([(1, "Riz"), (2, "Haricot")])
n_crops = length(crops)
successor = [[2], [1], [1,2]]
C = [[1],[2]] #1 = riz, 2 = haricot

#transitions[crop_before, crop_after, n_cultiv, n_jachere]

transitions = zeros(n_crops+1, n_crops+1, a_bar, l_bar)
transitions[n_crops+1, 1, 1, 1] = 72
transitions[n_crops+1, 1, 1, 2] = 120
transitions[n_crops+1, 2, 1, 1] = 54
transitions[n_crops+1, 2, 1, 2] = 90
transitions[1, 2, 2, 1] = 29
transitions[1, 2, 2, 2] = 65
transitions[2, 1, 2, 1] = 54
transitions[2, 1, 2, 2] = 90


demand = zeros(T,length(crops))
for time in 1:T
    demand[time, 2-time%2] = 400 * (1 + 2*(time%2))
end
display(demand)
n_states = (n_crops+1) * (n_crops+1) * l_bar

function arc_to_i_j(arc)
    t = arc[1]
    l = arc[2]
    a = arc[3]
    i = arc[4]
    j = arc[5]
    initial = (t,l,a,i)
    if time==0
        if i != 0 || l!=2 || a!=0
            return false
        end
    end
    if time==T
        if j!=0
            return false
        end
    end
    if j == 0
        return initial, (t+1, min(2, l+1*(i==0) + (i != 0)), 0, j)
    end
    if j == 1
        if a == a_bar || (t+1)%2 == 0 || i==1
            return false
        end
        return initial, (t+1, l, a+1, j)
    end
    if j == 2
        if a == a_bar || (t+1)%2 == 1 || i==2
            return false
        end
        return initial, (t+1, l, a+1, j)
    end
end

function predecesseurs(noeud)
    preds = []
    t = noeud[1]
    l = noeud[2]
    a = noeud[3]
    j = noeud[4]
    if t==0
        return false
    end
    for a_ante in 0:a_bar
        for l_ante in 1:l_bar
            for i in 0:n_crops
                res = arc_to_i_j([t-1, l_ante, a_ante, i, j])
                if res != false && res[1] == noeud
                    preds += [res[0]]
                end
            end
        end
    end
    return preds
end

m = Model(Gurobi.Optimizer)
@variable(m, 0 <= x[0:T, 1:l_bar, 1:a_bar, 0:n_crops, 0:n_crops] <= P)

for time in 0:T+1
    for l in 1:l_bar
        for a in 0:a_bar
            for crop_origin in 0:n_crops
                s_plus = 0
                s_minus = 0
                for predecesseur in predecesseurs((time, l, a, i))
                    s_minus += x[predecesseur[1], predecesseur[2], predecesseur[3], predecesseur[4], i]
                end
                for crop_destination in 0:n_crops
                    res = arc_to_i_j([time, l, a, crop_origin, crop_destination])
                    if res == false
                        @constraint(m, x[time, l, a, crop_origin, crop_destination] == 0)
                    else
                        origin, destination = res
                    end 
                    s_plus += x[time, l, a, crop_origin, crop_destination]
                end
                @constraint(m, s_minus - s_plus == 0)
            end
        end
    end
end

@constraint(m, [t in 1:T, crop in 1:n_crops], sum(x[t-1, a, l, i, crop]*transitions[i,crop+(n_crops+1)*(crop==0),a,l] for a in 0:a_bar, l in 1:l_bar, i in 0:n_crops) >= demand[t, crop])
 
@objective(m, Min, sum(x[0, l, a, 0, crop_dest] for a in 0:a_bar, l in 1:l_bar, crop_dest in 0:n_crops))

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