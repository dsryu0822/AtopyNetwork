include("../core/header.jl")

files = [(joinpath.(root, files) for (root, dirs, files) in walkdir("data"))...;]
filter!(x -> occursin(".txt", x), files)
feature_names = ["pareto", "max_central", "max_weight", "max_core", "max_ind", "max_pagerank", "min_cut", "min_dom", "diameter", "radius", "num_leaf", "coef_clustering", "coef_assortativity", "entropy", "mad", "runtime", ]

featureL = DataFrame(;
    ID = replace.(last.(split.(files, "\\")), ".txt" => ""),
    # y = occursin.("accept", files),
    y = .!(occursin.("reject", files) .|| occursin.("submitted", files) .|| occursin.("withdrawn", files)),
    year = parse.(Int64, getindex.(split.(files, "\\"), 2))
)

binary_check = ([240, 423, 500, 684, 856, 1091, 1570, 2257] .== [count(df.y) for df in groupby(featureL, :year)])
if all(binary_check)
    @info "binary check passed!"
else
    @error "binary check not passed: $(unique(featureL.year)[.!binary_check])"
end

featureR = DataFrame([zeros(length(files)) for _ in feature_names], feature_names)
feature = [featureL featureR]

# @warn "feature for year < 2022"
# feature = feature[feature.year .< 2022, :]
# feature = CSV.read("graph.csv", DataFrame)

drs = eachrow(feature)
errors = []
@showprogress @threads for idx in eachindex(drs)
try
    tic = now()
    dr = drs[idx]
    if dr.connected > 0 continue end
    file = files[idx]
    txt = read(file, String)
    crps = [txt |> StringDocument] |> Corpus
    remove_case!(crps)
    update_lexicon!(crps)

    spX = Int64.(CooMatrix(crps[1], window = 1, normalize = false).coom)
    for k in axes(spX, 1) spX[k, k] = 0 end
    G = Graph(spX)
    uG = UG(G)
    deg_G = degree(G)
    freq_G = [count(deg_G .== k) for k in sort(unique(deg_G))] / length(deg_G)


    # dr.y = ifelse(occursin("accept", file), 1, 0)
    dr.pareto = 1/fit(Pareto, deg_G).α
    dr.max_central = maximum(katz_centrality(G))
    dr.max_weight = maximum(spX)
    dr.max_core = maximum(k_core(G))
    dr.max_ind = length(independent_set(G, MaximalIndependentSet()))
    dr.max_pagerank = maximum(pagerank(G))
    dr.min_cut = only(min_cut(uG))
    dr.min_dom = length(dominating_set(G, MinimalDominatingSet()))
    dr.diameter = diameter(G)
    # dr.girth = girth(uG)
    dr.radius = Graphs.radius(G)
    # dr.connected = Graphs.is_connected(G)
    dr.num_leaf = count(deg_G .== 1)
    # dr.num_chromatic = chromatic_number(uG)
    dr.coef_clustering = global_clustering_coefficient(G)
    dr.coef_assortativity = assortativity(G)
    dr.entropy = Distributions.entropy(freq_G)
    dr.mad = SimpleGraphAlgorithms.mad(uG)
    dr.runtime = (now() - tic).value / 60000
    CSV.write("graph.csv", feature)
catch e
    push!(errors, idx)
end
end
feature = feature[Not(errors), :]
CSV.write("graph.csv", feature)


# @info "here are appendix:"
# 1/fit(Pareto, rand(1000) .^ (-3)).α
# @time betweenness_centrality(G);
# @time closeness_centrality(G);
# @time degree_centrality(G);
# @time katz_centrality(G);
# @time radiality_centrality(G);
# @time stress_centrality(G);
# Distributions.entropy([0.33, 0.33, 0.33])




# txt_ = read.(files, String)
# crps = Corpus(StringDocument.(txt_))
# remove_case!(crps)
# @time update_lexicon!(crps)

# @time spX = getproperty.(CooMatrix.(crps, window = 1, normalize = false), :coom)