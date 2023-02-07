using DataFrames
using CSV

include("./utils/constants.jl")
include("./models/modelPath.jl")

"""
include("./main.jl")
solveAllInstances()
"""

function runInstanceAndUpdateDataframe(currentResults::DataFrame, fileToRun::String, timeLimit::Float64, rowToReplace::Union{Int, Nothing}=nothing)::Bool
    
    result = pathSolve(fileToRun, timeLimit, true)
    if result == nothing
        println("NOT FEASIBLE!!")
        return false
    end

    optimal, solveTime, value, bound, gap = result
    # Modify dataframe
    if rowToReplace == nothing
        rowToReplace = findfirst(==(fileToRun), currentResults.instance)
        if rowToReplace == nothing
            println("Pushing new row to results dataframe")
            push!(currentResults, [fileToRun optimal solveTime value bound gap])
            return true
        else
            currentRow= currentResults[rowToReplace,:]
            if value < currentRow.value || (solveTime > currentRow.time && value < currentRow.value + 1e-6)
                println("Improved value for " * fileToRun)
                currentResults[rowToReplace,:] = [fileToRun optimal solveTime value bound gap]
                return true
            end
        end
    else
        currentRow = currentResults[rowToReplace,:]
        if value < currentRow.value || (solveTime > currentRow.time && value < currentRow.value + 1e-6)
            println("Improved value for " * fileToRun)
            currentResults[rowToReplace,:] = [fileToRun optimal solveTime value bound gap]
            return true
        end
    end
    return false
end


function solveAllInstances(resultFile::String=RESULTS_FILE, timeLimit::Float64=-1.)::Nothing
    # Loading
    filePath =RESULTS_DIR_PATH * "\\" * resultFile * ".csv"
    # Get unoptimal instance
    if !isfile(filePath)
        currentResults = DataFrame(instance=String[], optimal=Bool[], time=Float64[], value=Float64[], bound=Float64[], gap=Float64[])
    else
        currentResults = DataFrame(CSV.File(filePath))
    end

    # Run
    for fileToRun in DATA_FILES
        updatedDf = runInstanceAndUpdateDataframe(currentResults, fileToRun, timeLimit)
        if updatedDf
            CSV.write(filePath, currentResults, delim=";")
        end
    end
    return 
end
