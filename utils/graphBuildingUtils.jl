using Graphs

function buildGraphFromInstance(n::Int64, A_1::Matrix{Int64}, A_2::Matrix{Int64})::DiGraph
    graph = DiGraph(n)
    for (i,j,_) in eachrow(A_1) 
        add_edge!(graph, i, j)
    end

    for (i,j,_) in eachrow(A_2)
        add_edge!(graph, i, j)
    end

    return graph
end