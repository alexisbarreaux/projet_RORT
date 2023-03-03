using JuMP
using CPLEX

include("../utils/constants.jl")
include("../utils/graphBuildingUtils.jl")
include("../utils/jsonUtils.jl")
include("../models/modelPath.jl")

"""
include("./heuristics/heuristic3.jl")
heurConvexSolve("taxe_grille_2x3.txt")
"""

function heurConvexSolve(inputFile::String)
    T_val, relaxedSolveTime = pathSolve(inputFile, relaxed = true, boundMode=2, convex = true)
    isOptimal, fixedSolveTime, value, bound, gap= pathSolve(inputFile, T_val = T_val, boundMode=2)
    return isOptimal, fixedSolveTime + relaxedSolveTime, value, bound, gap
end