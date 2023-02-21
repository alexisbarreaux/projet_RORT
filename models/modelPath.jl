using JuMP
using CPLEX

include("../utils/constants.jl")
include("../utils/graphBuildingUtils.jl")
include("../utils/jsonUtils.jl")

"""
include("./models/modelPath.jl")
pathSolve("taxe_grille_2x3.txt")
pathSolve("dummy_graph.txt")

# Expected result on dummy_graph is 3
"""

function pathSolve(inputFile::String; timeLimit::Float64= -1., silent::Bool=true, boundMode::Int64=2, relaxed::Bool=false, T_val::Any=nothing)::Any
    """
    """
    if relaxed
        println("Solving relaxation")
    elseif T_val != nothing
        println("Solving with fixed T")
    else
        printtln("Solving ", inputFile)
    end
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

    if boundMode == 1
        M = sum(A_2[i,3] for i in 1:numberOfArcs_A2)
    elseif boundMode == 2
        A2OnlyGraph = buildA2Graph(n, A_2)
        
        M = Dict{Tuple{Int64, Int64}, Float64}()
        weights = Matrix{Float64}(undef, n, n)
        for (i, j, c_a) in eachrow(A_2)
            weights[i,j] = c_a
        end

        # Measuring M_a_k
        for (s_k, t_k) in eachrow(K)
            ds = dijkstra_shortest_paths(A2OnlyGraph, s_k, weights)
            M[(s_k, t_k)] = ds.dists[t_k]
        end

        # Measuring
        M_a =  Dict{Tuple{Int64, Int64}, Float64}()
        for (i, j, c_a) in eachrow(A_1)
            M_a[(i,j)] = maximum([M[(s_k, t_k)] - c_a for (s_k, t_k) in eachrow(K)])
        end
    else
        M = 100
    end

    ### Variables
    # Fix T if needed
    if T_val != nothing
        T = T_val
    else
        @variable(model, T[i in 1:n, j in 1:n] >= 0.0)
    end
    if !relaxed
        @variable(model, y[i in 1:n, j in 1:n, k in 1:numberOfCommodities], Bin)
    else
        @variable(model, 0. <= y[i in 1:n, j in 1:n, k in 1:numberOfCommodities] <= 1.)
    end
    @variable(model, z[i in 1:n, j in 1:n, k in 1:numberOfCommodities] >= 0.0)
    # Variables from second level
    @variable(model, alpha[i in 1:n, k in 1:numberOfCommodities])

    # Objective : sum on a in A_1 and k in K of the z_a * n_k
    @objective(model, Max, sum(z[i, j, k] * n_k[k] for (i, j, _) in eachrow(A_1) for k in 1:numberOfCommodities))

    ### Constraints
    # Linearisation constraints
    @constraint(model, [(i, j, _) in eachrow(A_1), k in 1:numberOfCommodities], z[i, j, k] <= T[i, j])
    if boundMode != 2
        @constraint(model, [(i, j, _) in eachrow(A_1), k in 1:numberOfCommodities], z[i, j, k] <= y[i, j, k] * M)
        @constraint(model, [(i, j, _) in eachrow(A_1), k in 1:numberOfCommodities], T[i, j] - z[i, j, k] <= M * (1 - y[i, j, k]))
    else
        @constraint(model, [(i, j, c_a) in eachrow(A_1), k in 1:numberOfCommodities], z[i, j, k] <= y[i, j, k] * max(M[(K[k,1], K[k,2])] - c_a, 0))
        @constraint(model, [(i, j, c_a) in eachrow(A_1), k in 1:numberOfCommodities], T[i, j] - z[i, j, k] <= max(M_a[(i,j)], 0)  * (1 - y[i, j, k]))
    end

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
        @constraint(model, alpha[source, k] + alpha[dest, k] >= sum(c_a * y[i, j, k] + z[i, j, k] for (i, j, c_a) in eachrow(A_1)) + sum(c_a * y[i, j, k] for (i, j, c_a) in eachrow(A_2)))

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

    if !silent
        y_val = JuMP.value.(y)
        T_val = JuMP.value.(T)
        z_val = JuMP.value.(z)
        for (i, j, _) in eachrow(A_1)
            if T_val[i, j] > 0.0
                println("Edge " * string(i) * " -> " * string(j) * " toll is " * string(T_val[i, j]))
                for k in 1:numberOfCommodities
                    if y_val[i, j, k] > 0.0
                        println("Commodity " * string(k) * " uses edge. z is " * string(z_val[i, j, k]))
                    end
                end
            end
        end
    end

    
    if feasibleSolutionFound
        solveTime = round(JuMP.solve_time(model), digits=5)
        if relaxed
            return JuMP.value.(T), solveTime
        else
            value = round(JuMP.objective_value(model), digits=5)
            gap = JuMP.relative_gap(model)
            bound = JuMP.objective_bound(model)
            return isOptimal, solveTime, value, bound, gap
        end
    else
        println("Problem is not feasible !!!")
        return
    end
end