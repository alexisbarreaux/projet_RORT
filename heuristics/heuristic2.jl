using JuMP
using CPLEX

include("../utils/constants.jl")
include("../utils/graphBuildingUtils.jl")
include("../utils/jsonUtils.jl")
include("../models/modelPath.jl")
include("../models/evalModel.jl")
include("../models/evalModel2.jl")

"""
include("./heuristics/heuristic2.jl")
heurSolveLocal2("taxe_grille_5x6.txt")
heurSolveLocal("taxe_grille_5x6.txt")
heurSolve("taxe_grille_2x3.txt")
"""

function heurSolve(inputFile::String)
    T_val, relaxedSolveTime = pathSolve(inputFile, relaxed = true, boundMode=2)
    isOptimal, fixedSolveTime, value, bound, gap= evalSolve(inputFile, T_val = T_val)
    return isOptimal, fixedSolveTime + relaxedSolveTime, value, bound, gap
end

function heurSolveLocal(inputFile::String)
    start = time()
    T_val, relaxedSolveTime, M_a, M_a_k, numberOfCommodities = pathSolve(inputFile, relaxed = true, boundMode=3, heur2=true)
    # First solve
    isOptimal, _, bestValue, _, _= evalSolve(inputFile, T_val = T_val)

    # Try to add more toll on arc with biggest gap between toll and bound
    # For each arc possible to tax
    A_1_arcs = keys(M_a)

    for (i,j) in A_1_arcs
        # desc sort steps
        steps = sort(vcat(unique([M_a_k[(i,j,k)] for k in 1:numberOfCommodities]), 0), rev=true)
        for step in steps
            if T_val[i,j] <= step
                continue
            else
                oldValue = T_val[i,j]
                T_val[i,j] = step
                _, _, newValue, _, _= evalSolve(inputFile, T_val = T_val)
                if bestValue <= newValue
                    bestValue = newValue
                    continue
                else
                    T_val[i,j] = oldValue
                    break
                end
            end
        end

    end

    return isOptimal, time() - start, bestValue, 0., 0.
end

function heurSolveLocal2(inputFile::String)
    start = time()
    T_val, relaxedSolveTime, M_a, M_a_k, numberOfCommodities = pathSolve(inputFile, relaxed = true, boundMode=3, heur2=true)
    # First solve
    model, dictA_1, dictA_2, n_k, numberOfCommodities = initModel(inputFile, T_val)

    isOptimal, _, bestValue, _, _ = evaluateModel(model)
    # Try to add more toll on arc with biggest gap between toll and bound
    # For each arc possible to tax
    A_1_arcs = keys(M_a)

    for (i,j) in A_1_arcs
        # desc sort steps
        steps = sort(vcat(unique([M_a_k[(i,j,k)] for k in 1:numberOfCommodities]), 0), rev=true)
        for step in steps
            if T_val[i,j] <= step
                continue
            else
                oldValue = T_val[i,j]
                T_val[i,j] = step
                modifyModel(model, i, j, T_val[i,j], dictA_1, dictA_2, n_k, numberOfCommodities)
                _, _, newValue, _, _= evaluateModel(model)
                if bestValue <= newValue
                    bestValue = newValue
                    continue
                else
                    T_val[i,j] = oldValue
                    modifyModel(model, i, j, T_val[i,j], dictA_1, dictA_2, n_k, numberOfCommodities)
                    break
                end
            end
        end

    end

    return isOptimal, time() - start, bestValue, 0., 0.
end