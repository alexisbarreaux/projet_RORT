using JuMP
using CPLEX

include("../utils/constants.jl")
include("../utils/graphBuildingUtils.jl")

"""
include("./models/evalModel2.jl")
evalSolve2("taxe_grille_2x3.txt")
"""

function initModel(inputFile::String, T_val::Any)::Any
    # Directly load data file
    include(DATA_DIR_PATH * "\\" * inputFile)
    graph = buildGraphFromInstance(n, A_1, A_2)

    numberOfArcs_A1 = size(A_1)[1]
    numberOfArcs_A2 = size(A_2)[1]
    totalArcs = numberOfArcs_A1 + numberOfArcs_A2
    numberOfCommodities = size(K)[1]

    # Creating the model
    model = Model(CPLEX.Optimizer)
    set_silent(model)

    ### Variables
    @variable(model, 0. <= y[i in 1:n, j in 1:n, k in 1:numberOfCommodities] <= 1.)

    # Variables from second level
    @variable(model, alpha[i in 1:n, k in 1:numberOfCommodities])

    # Objective : sum on a in A_1 and k in K of the z_a * n_k
    @objective(model, Max, sum(T_val[i, j]*y[i, j, k] * n_k[k] for (i, j, _) in eachrow(A_1) for k in 1:numberOfCommodities))

    ### Constraints
    
    
    # Constraints from the sub problems
    for k in 1:numberOfCommodities
        (source, dest) = K[k, :]
        ## Primal constraints
        for i in 1:n
            if i == source
                @constraint(model, sum(y[i, succ, k] for succ in outneighbors(graph, i)) == 1)
            elseif i == dest
                @constraint(model, sum(y[pred, i, k] for pred in inneighbors(graph, i)) == 1)
            else
                @constraint(model, (sum(y[pred, i, k] for pred in inneighbors(graph, i)) - sum(y[i, succ, k] for succ in outneighbors(graph, i))) == 0)
            end
        end

        ## Dual constraints
        # Dual objective
        @constraint(model, alpha[source, k] + alpha[dest, k] >= sum( (c_a + T_val[i, j]) * y[i, j, k] for (i, j, c_a) in eachrow(A_1)) + sum(c_a * y[i, j, k]
         for (i, j, c_a) in eachrow(A_2)), base_name="dual_obj_"*string(k))

        for (i, j, c_a) in eachrow(A_1)
            # i = s
            if i == source
                # (s,t)
                if j == dest
                    @constraint(model, alpha[source, k] + alpha[dest, k] <= c_a + T_val[i, j], base_name=string(i)*"_"*string(j)*"_"*string(k)* "s,t")
                    # (s,j)
                else
                    @constraint(model, alpha[source, k] + alpha[j, k] <= c_a + T_val[i, j], base_name= string(i)*"_"*string(j)*"_"*string(k)* "s,j")
                    
                end
                # i = t
            elseif i == dest
                # (t, s)
                if j == source
                    continue
                    # (t,j)
                else
                    @constraint(model, alpha[j, k] <= c_a + T_val[i, j], base_name=string(i)*"_"*string(j)*"_"*string(k)* "t,j")
                    
                end
            else
                # (i,s)
                if j == source
                    @constraint(model, -alpha[i, k] <= c_a + T_val[i, j], base_name=string(i)*"_"*string(j)*"_"*string(k)* "i,s")
                    
                    # (i,t)
                elseif j == dest
                    @constraint(model,  alpha[dest, k] - alpha[i, k] <= c_a + T_val[i, j], base_name=string(i)*"_"*string(j)*"_"*string(k)* "i,t")
                    
                    # (i,j)
                else
                    @constraint(model, alpha[j, k] - alpha[i, k] <= c_a + T_val[i, j], base_name=string(i)*"_"*string(j)*"_"*string(k)* "i,j")
                    
                end
            end
        end

        for (i, j, c_a) in eachrow(A_2)
            # i = s
            if i == source
                # (s,t)
                if j == dest
                    @constraint(model, alpha[i, k] + alpha[j, k] <= c_a)
                    # (s,j)
                else
                    @constraint(model, alpha[i, k] + alpha[j, k] <= c_a)
                end
                # i = t
            elseif i == dest
                # (t, s)
                if j == source
                    continue
                    # (t,j)
                else
                    @constraint(model, alpha[j, k] <= c_a)
                end
            else
                # (i,s)
                if j == source
                    @constraint(model, -alpha[i, k] <= c_a)
                    # (i,t)
                elseif j == dest
                    @constraint(model, alpha[j, k] - alpha[i, k] <= c_a)
                    # (i,j)
                else
                    @constraint(model, alpha[j, k] - alpha[i, k] <= c_a)
                end
            end
        end
    end
    dictA_1=Dict{Tuple{Int64, Int64}, Int64}()
    dictA_2=Dict{Tuple{Int64, Int64}, Int64}()
    for (i,j,c_a) in eachrow(A_1)
        dictA_1[(i,j)] = c_a
    end
    for (i,j,c_a) in eachrow(A_2)
        dictA_2[(i,j)] = c_a
    end
    return model, dictA_1, dictA_2, n_k, numberOfCommodities
end

function modifyModel(model::Any, i::Int64, j::Int64, new_T_val::Any, dictA_1::Dict, dictA_2::Dict, n_k, numberOfCommodities)
    y = model[:y]

    for k in 1:numberOfCommodities
        set_objective_coefficient(model, y[i,j,k], new_T_val * n_k[k] )

        if haskey(dictA_1, (i,j))
            # Normalized coeff are sometimes the opposite of ours.
            c  =constraint_by_name(model, "dual_obj_"*string(k))

            if normalized_coefficient(c, y[i,j,k]) < 0
                newCoeff = - (new_T_val + dictA_1[(i,j)])
            else
                newCoeff = new_T_val + dictA_1[(i,j)]
            end
            set_normalized_coefficient(c, y[i,j,k], newCoeff)
        end

        constraint = constraint_by_name(model, string(i)*"_"*string(j)*"_"*string(k)* "s,t")
        if constraint != nothing
            if normalized_rhs(constraint) < 0
                newCoeff = - (new_T_val + dictA_1[(i,j)])
            else
                newCoeff = new_T_val + dictA_1[(i,j)]
            end
            set_normalized_rhs(constraint, newCoeff)
        end
        constraint = constraint_by_name(model, string(i)*"_"*string(j)*"_"*string(k)* "s,j")
        if constraint != nothing
            if normalized_rhs(constraint) < 0
                newCoeff = - (new_T_val + dictA_1[(i,j)])
            else
                newCoeff = new_T_val + dictA_1[(i,j)]
            end
            set_normalized_rhs(constraint, newCoeff)
        end
        constraint = constraint_by_name(model, string(i)*"_"*string(j)*"_"*string(k)* "t,j")
        if constraint != nothing
            if normalized_rhs(constraint) < 0
                newCoeff = - (new_T_val + dictA_1[(i,j)])
            else
                newCoeff = new_T_val + dictA_1[(i,j)]
            end
            set_normalized_rhs(constraint, newCoeff)
        end
        constraint = constraint_by_name(model, string(i)*"_"*string(j)*"_"*string(k)* "i,s")
        if constraint != nothing
            if normalized_rhs(constraint) < 0
                newCoeff = - (new_T_val + dictA_1[(i,j)])
            else
                newCoeff = new_T_val + dictA_1[(i,j)]
            end
            set_normalized_rhs(constraint, newCoeff)
        end
        constraint = constraint_by_name(model, string(i)*"_"*string(j)*"_"*string(k)* "i,t")
        if constraint != nothing
            if normalized_rhs(constraint) < 0
                newCoeff = - (new_T_val + dictA_1[(i,j)])
            else
                newCoeff = new_T_val + dictA_1[(i,j)]
            end
            set_normalized_rhs(constraint, newCoeff)
        end
        constraint = constraint_by_name(model, string(i)*"_"*string(j)*"_"*string(k)* "i,j")
        if constraint != nothing
            if normalized_rhs(constraint) < 0
                newCoeff = - (new_T_val + dictA_1[(i,j)])
            else
                newCoeff = new_T_val + dictA_1[(i,j)]
            end
            set_normalized_rhs(constraint, newCoeff)
        end
    end
    return
end


function evaluateModel(model::Any)
    optimize!(model)
    feasibleSolutionFound = primal_status(model) == MOI.FEASIBLE_POINT
    isOptimal = termination_status(model) == MOI.OPTIMAL
    
    if feasibleSolutionFound
        solveTime = round(JuMP.solve_time(model), digits=5)
        value = round(JuMP.objective_value(model), digits=5)
        return isOptimal, solveTime, value, 0., 0.
    else
        println("Problem is not feasible !!!")
        return
    end
end
