using JuMP
using CPLEX

include("../utils/constants.jl")
include("../utils/graphBuildingUtils.jl")
include("../utils/jsonUtils.jl")
include("../models/modelPath.jl")

"""
include("./heuristics/heuristic2.jl")
heurSolve("taxe_grille_2x3.txt")
"""

function heurSolve(inputFile::String)
    T_val, relaxedSolveTime = pathSolve(inputFile, relaxed = true, boundMode=2)
    isOptimal, fixedSolveTime, value, bound, gap= pathSolve(inputFile, T_val = T_val, boundMode=2)
    return isOptimal, fixedSolveTime + relaxedSolveTime, value, bound, gap
end