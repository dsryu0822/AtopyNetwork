using ProgressMeter
packages = [:JLD2, :SparseArrays, :DataFrames, :Plots, :Random, :TextAnalysis,
 :Snowball, :CSV, :Graphs, :GraphRecipes, :LinearAlgebra, :NetworkLayout,
 :StatsBase, :LaTeXStrings]
@showprogress for package in packages
    @eval using $(package)
end
mm = Plots.mm

struct WordNetwork
    adj::AbstractSparseMatrix
    nodes::Vector{String}
end
function Base.show(io::IO, WN::WordNetwork)
    print(io, "$(length(WN.nodes))-node WordNetwork: ")
    show(io, WN.nodes)
end

function neighbor(WN::WordNetwork, id_query::Integer)
    if !(1 ≤ id_query ≤ length(WN.nodes))
        error(BoundsError, ": node with id '$id_query' is not in the network")
    end

    id = findall(.!iszero.(WN.adj[id_query, :]))
    node = WN.nodes[id]
    weight = WN.adj[id_query, id]

    df = DataFrame(; id, node, weight)
    sort!(df, :weight, rev = true)
    return df
end
function neighbor(WN::WordNetwork, query::String)
    id_query = findfirst(WN.nodes .== query)
    if id_query |> isnothing
        error(KeyError, ": node with name '$query' is not in the network")
    end
    return neighbor(WN, id_query)
end

function subgraph(WN::WordNetwork, id, depth = 5, distance = 1)
    if id isa String
        id = findfirst(WN.nodes .== id)
    end

    df  = DataFrame(depth = 0, src = 0, id = id, weight = 0)
    for i in 1:distance
        for j in df[df.depth .== (i-1), :id]
            if j in df.src continue end
            df = [df; [
                DataFrame(depth = fill(i, depth), src = j) neighbor(WN, j)[1:depth, [:id, :weight]]
                ]]
        end
    end
    rename!(df, :id => :dst)
    df = df[df.src .!= df.dst, :]
    # unique!(df, :dst)
    return df
end

function plot_subgraph(WN::WordNetwork, id, depth = 5, distance = 2)
    distance = distance + 2
    SG = subgraph(WN, id, depth, distance)
    _id = setdiff(unique([SG.src; SG.dst]), 0)
    _node = WN.nodes[_id]
    sn = length(_id)
    
    dict_id = Dict(_id .=> 1:sn)
    
    θ = shuffle(LinRange(0, 2π, sn) + 0.01randn(sn))
    r = SG.depth[[findfirst(SG.dst .== id) for id in _id]]
    x_ = r .* cos.(θ)
    y_ = r .* sin.(θ)
    
    SG.src[1] = first(SG.dst)
    srcx = x_[replace(SG.src, dict_id...)]
    srcy = y_[replace(SG.src, dict_id...)]
    dstx = x_[replace(SG.dst, dict_id...)]
    dsty = y_[replace(SG.dst, dict_id...)]
    
    p1 = plot()
    for mat = eachrow(stack([srcx, dstx, srcy, dsty]))
        plot!(p1, mat[1:2], mat[3:4], color = :black, lw = 0.5, alpha = 0.5)
    end
    p1 = scatter!(x_, y_, text = _node,
        color = :white, msw = 0, ms = 20, size = (900, 900), legend = :none, alpha = 0.8,
        lims = (-distance+1, distance-1), framestyle=:none)

    return p1
end


function logscatter(degree; scale = :identity, xscale = :identity, yscale = :identity, kargs...)
    if scale == :log10 || xscale == :log10
        xtick = exp10.(0:.1:ceil(log10(maximum(degree))))
    else
        xtick = eps():50:ceil(maximum(degree))
    end
    ytick = zeros(Int64, length(xtick))
    for dg = degree
        idx = findfirst(dg .≤ xtick)
        if !isnothing(idx)
            ytick[findfirst(dg .≤ xtick)] += 1
        end
    end
    xtick, ytick = xtick[.!iszero.(ytick)], ytick[.!iszero.(ytick)]
    if xscale == :log10
        return scatter(xtick, ytick; xscale, kargs...)
    elseif yscale == :log10
        return scatter(xtick, ytick; yscale, kargs...)
    else
        return scatter(xtick, ytick; scale, kargs...)
    end
end