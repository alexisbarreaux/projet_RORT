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

function solveAll(resultFile::String="results.json")
    results =  Dict{String, Float64}()
    for file in DATA_FILES
        results[file] = pathSolve(file)
    end
    filePath =RESULTS_DIR_PATH * "\\" * resultFile
    jsonDropToFile(filePath, results)
end

function pathSolve(inputFile::String, silent::Bool=true)::Any
    """
    """
    println("Solving ", inputFile)
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

    M = 100 # TODO improve the big M

    ### Variables
    @variable(model, T[i in 1:n, j in 1:n] >= 0)
    @variable(model, y[i in 1:n, j in 1:n, k in 1:numberOfCommodities], Bin)
    @variable(model, z[i in 1:n, j in 1:n, k in 1:numberOfCommodities] >= 0)
    # Variables from second level
    @variable(model, alpha[i in 1:n, k in 1:numberOfCommodities])
    
    # Objective : sum on a in A_1 and k in K of the z_a * n_k
    @objective(model, Max, sum(z[i,j,k] * n_k[k] for (i,j,_) in eachrow(A_1) for k in 1:numberOfCommodities))

    ### Constraints
    # Linearisation constraints
    @constraint(model, [(i,j,_) in eachrow(A_1), k in 1:numberOfCommodities], z[i,j,k] <= T[i,j])
    @constraint(model, [(i,j,_) in eachrow(A_1), k in 1:numberOfCommodities], z[i,j,k] <= y[i,j,k]*M)

    # Constraints from the sub problems
    for k in 1:numberOfCommodities
        (source, dest) = K[k,:]
        ## Primal constraints
        for i in 1:n
            if i == source
                @constraint(model, sum(y[i,succ,k] for succ in outneighbors(graph, i)) == 1)
            elseif i == dest
                @constraint(model, sum(y[pred,i,k] for pred in inneighbors(graph, i)) == 1)
            else
                @constraint(model, (sum(y[pred,i,k] for pred in inneighbors(graph, i)) - sum(y[i,succ,k] for succ in outneighbors(graph, i))) == 0)
            end
        end
        ## Dual constraints
        # Dual objective
        @constraint(model, alpha[source,k] + alpha[dest,k] >= sum(c_a*y[i,j,k] + z[i,j,k] for (i,j,c_a) in eachrow(A_1)) + sum(c_a*y[i,j,k] for (i,j,c_a) in eachrow(A_2)))

        for (i,j,c_a) in eachrow(A_1)
            if i == source
                @constraint(model, alpha[i,k] + alpha[j,k] <= c_a + T[i,j])
            elseif i == dest
                if j == source
                    continue
                else 
                    @constraint(model, alpha[j,k] <= c_a + T[i,j])
                end
            else
                if j == source
                    @constraint(model, -alpha[i,k] <= c_a + T[i,j])
                else
                    @constraint(model, alpha[j,k] -alpha[i,k] <= c_a + T[i,j])
                end
            end
        end
    
        for (i,j,c_a) in eachrow(A_2)
            if i == source
                @constraint(model, alpha[i,k] + alpha[j,k] <= c_a)
            elseif i == dest
                if j == source
                    continue
                else 
                    @constraint(model, alpha[j,k] <= c_a)
                end
            else
                if j == source
                    @constraint(model, -alpha[i,k] <= c_a)
                else
                    @constraint(model, alpha[j,k] - alpha[i,k] <= c_a)
                end
            end
        end
    end


    optimize!(model)
    feasibleSolutionFound = primal_status(model) == MOI.FEASIBLE_POINT
    isOptimal = termination_status(model) == MOI.OPTIMAL
    if feasibleSolutionFound
        return JuMP.objective_value(model)
    else
        println("Problem is not feasible !!!")
        return
    end
end