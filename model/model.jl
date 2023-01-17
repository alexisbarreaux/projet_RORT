using JuMP
using CPLEX


function exactSolve(inputFile::String, showResult::Bool= false, silent::Bool=true)::Any
    """
    """
    println("Solving ", inputFile, " in exact mode.")
    # Directly load data file
    include(DATA_DIR_PATH * "\\" * inputFile)

    l = computeDistances(coordinates)

    # Creating the model
    model = Model(CPLEX.Optimizer)
    if silent
        set_silent(model)
    end

    # Variables
    @variable(model, x[i in 1:n, j in 1:n], Bin)
    @variable(model, y[i in 1:n, k in 1:K], Bin)

    # Objective
    # We only sum on i < j as per how the subject specified it on an example
    @objective(model, Min, sum(x[i,j] * l[i,j] for i in 1:n for j in 1:n))
    
    # Constraints
    # Weights of the parts
    @constraint(model, [k in 1:K], sum(w_v[i]*y[i,k] for i in 1:n) <= B)
    # Triangular inequalities between x and y
    @constraint(model, [k in 1:K, i in 1:n, j in i:n], y[i,k] + y[j,k] <= x[i,j] + 1)
    # Each node is in a part
    @constraint(model, [i in 1:n],  sum(y[i,k] for k in 1:K) == 1)

    # Solve
    start = time()
    optimize!(model)
    optimize_time = time() - start

    ### Display the solution
    feasibleSolutionFound = primal_status(model) == MOI.FEASIBLE_POINT
    isOptimal = termination_status(model) == MOI.OPTIMAL
    if feasibleSolutionFound
    # Récupération des valeurs d’une variable
        result = JuMP.value.(y)
        value = JuMP.objective_value(model)
        
        if showResult
            println("Success, nodes : " * string(JuMP.node_count(model))* ", Time : "* string(round(optimize_time, digits= 5)) * " Value : " * string(round(value, digits=4)))
            createdParts = Dict{Int, Array}(k => [] for k in 1:K)
            for i in 1:n
                for k in 1:K
                    if result[i,k] == 1
                        createdParts[k] = vcat(createdParts[k],[i])
                        break
                    end
                end
            end
            println("Found parts are : ", createdParts)
        end
        return value
    else
        println("Not feasible!!")
        return
    end

end