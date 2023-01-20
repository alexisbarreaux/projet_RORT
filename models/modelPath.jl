using JuMP
using CPLEX

include("../utils/constants.jl")

"""
include("./models/modelPath.jl")
pathSolve("taxe_grille_2x3.txt")
pathSolve("dummy_graph.txt")

# Expected result on dummy_graph is 3
"""


function pathSolve(inputFile::String, silent::Bool=true)::Any
    """
    """
    println("Solving ", inputFile)
    # Directly load data file
    include(DATA_DIR_PATH * "\\" * inputFile)

    numberOfArcs_A1 = size(A_1)[1]
    numberOfArcs_A2 = size(A_2)[1]
    totalArcs = numberOfArcs_A1 + numberOfArcs_A2
    numberOfCommodities = size(K)[1]

    function kThSubProblem(k::Int64)::Any
        # Creating the model
        model = Model(CPLEX.Optimizer)
        if silent
            set_silent(model)
        end

        (source, dest) = K[k,:]
        M = 100 # TODO improve the big M

        # Variables
        @variable(model, T[i in 1:numberOfArcs_A1] >= 0)
        @variable(model, y[i in 1:totalArcs], Bin)
        @variable(model, z[i in 1:numberOfArcs_A1] >= 0)
        # Variables from second level
        @variable(model, alpha[i in 1:n])
    
        # Objective : sum on a in A_1 and k in K of the z_a * n_k
        @objective(model, Max, sum(z[i] * n_k[k] for i in 1:numberOfArcs_A1))
    
        # Constraints
        # Linearisation constraints
        @constraint(model, [i in 1:numberOfArcs_A1], z[i] <= T[i])
        @constraint(model, [i in 1:numberOfArcs_A1], z[i] <= y[i]*M)
        # Second level induced constraints        
        @constraint(model, alpha[source] + alpha[dest] >= sum(A_1[i,3]*y[i] + z[i] for i in 1:numberOfArcs_A1) + sum(A_2[i,3]*y[i] for i in 1:numberOfArcs_A2))
    
        for a in 1:numberOfArcs_A1
            (i,j,c) = A_1[a,:]
            if i == source
                @constraint(model, alpha[i] + alpha[j] <= c + T[a])
            elseif i == dest
                if j == source
                    continue
                else 
                    @constraint(model, alpha[j] <= c + T[a])
                end
            else
                if j == source
                    @constraint(model, -alpha[i] <= c + T[a])
                else
                    @constraint(model, alpha[j] -alpha[i] <= c + T[a])
                end
            end
        end

        for a in 1:numberOfArcs_A2
            (i,j,c) = A_2[a,:]
            if i == source
                @constraint(model, alpha[i] + alpha[j] <= c)
            elseif i == dest
                if j == source
                    continue
                else 
                    @constraint(model, alpha[j] <= c)
                end
            else
                if j == source
                    @constraint(model, -alpha[i] <= c)
                else
                    @constraint(model, alpha[j] - alpha[i] <= c)
                end
            end
        end

        optimize!(model)
        feasibleSolutionFound = primal_status(model) == MOI.FEASIBLE_POINT
        isOptimal = termination_status(model) == MOI.OPTIMAL
        if feasibleSolutionFound
            return JuMP.objective_value(model)
        else
            println("Sub problem " * string(k) * " is not feasible !!!")
            return
        end

    end

    result = 0
    
    for k in 1:numberOfCommodities
        # n_k[k] already included in subproblem
        result += kThSubProblem(k)
    end
    return result
end