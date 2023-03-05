using JuMP
using CPLEX

include("../utils/constants.jl")
include("../utils/graphBuildingUtils.jl")
include("../utils/jsonUtils.jl")

"""
include("./models/evalModel.jl")
evalSolve("taxe_grille_2x3.txt")

# Expected result on dummy_graph is 3
"""

function evalSolve(inputFile::String; timeLimit::Float64= -1., silent::Bool=true, T_val::Any=nothing)::Any
    """
    """
    println("Solving with fixed T")

    # Directly load data file
    include(DATA_DIR_PATH * "\\" * inputFile)
    graph = buildGraphFromInstance(n, A_1, A_2)
    
    numberOfArcs_A1 = size(A_1)[1]
    numberOfArcs_A2 = size(A_2)[1]
    totalArcs = numberOfArcs_A1 + numberOfArcs_A2
    numberOfCommodities = size(K)[1]

    # Creating the model
    model = Model(CPLEX.Optimizer)
    if silent
        set_silent(model)
    end
    if timeLimit >= 0
        set_time_limit_sec(model, timeLimit)
    end

    ### Variables
    # Fix T
    T = T_val
    @variable(model, y[i in 1:n, j in 1:n, k in 1:numberOfCommodities], Bin)

    # Variables from second level
    @variable(model, alpha[i in 1:n, k in 1:numberOfCommodities])

    # Objective : sum on a in A_1 and k in K of the z_a * n_k
    @objective(model, Max, sum(T[i, j]*y[i, j, k] * n_k[k] for (i, j, _) in eachrow(A_1) for k in 1:numberOfCommodities))

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
        @constraint(model, alpha[source, k] + alpha[dest, k] >= sum( (c_a + T[i, j]) * y[i, j, k] for (i, j, c_a) in eachrow(A_1)) + sum(c_a * y[i, j, k] for (i, j, c_a) in eachrow(A_2)))

        for (i, j, c_a) in eachrow(A_1)
            # i = s
            if i == source
                # (s,t)
                if j == dest
                    @constraint(model, alpha[source, k] + alpha[dest, k] <= c_a + T[i, j])
                    # (s,j)
                else
                    @constraint(model, alpha[source, k] + alpha[j, k] <= c_a + T[i, j])
                end
                # i = t
            elseif i == dest
                # (t, s)
                if j == source
                    continue
                    # (t,j)
                else
                    @constraint(model, alpha[j, k] <= c_a + T[i, j])
                end
            else
                # (i,s)
                if j == source
                    @constraint(model, -alpha[i, k] <= c_a + T[i, j])
                    # (i,t)
                elseif j == dest
                    @constraint(model, alpha[dest, k] - alpha[i, k] <= c_a + T[i, j])
                    # (i,j)
                else
                    @constraint(model, alpha[j, k] - alpha[i, k] <= c_a + T[i, j])
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


    optimize!(model)
    feasibleSolutionFound = primal_status(model) == MOI.FEASIBLE_POINT
    isOptimal = termination_status(model) == MOI.OPTIMAL
    
    if feasibleSolutionFound
        solveTime = round(JuMP.solve_time(model), digits=5)
        value = round(JuMP.objective_value(model), digits=5)
        gap = JuMP.relative_gap(model)
        bound = JuMP.objective_bound(model)
        return isOptimal, solveTime, value, bound, gap
    else
        println("Problem is not feasible !!!")
        return
    end
end