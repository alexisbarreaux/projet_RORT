import JSON

function jsonDropToFile(filePath::String,data::Any)::Nothing
    jsonData = JSON.json(data)
    open(filePath, "w") do file
        write(file, jsonData)
    end
    return
end

function jsonLoadDictFromFile(filePath::String)::Dict
    open(filePath, "r") do file
        parsedDict = JSON.parsefile(filePath; dicttype=Dict{String, Float64}, inttype=Int64, use_mmap=true)
        return parsedDict
    end
end
