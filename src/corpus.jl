using ProgressMeter
packages = [:JLD2, :SparseArrays, :DataFrames, :Plots, :Random, :TextAnalysis,
 :Snowball, :CSV, :Graphs, :GraphRecipes, :LinearAlgebra, :NetworkLayout]
@showprogress for package in packages
    @eval using $(package)
end
mm = Plots.mm

cd("G:/ElderlyAD/")

gowords = CSV.read("gowords.csv", DataFrame)[:, 1]
stopwords = CSV.read("./clinical-concepts/clinical-stopwords.txt", DataFrame)[:, 1]
append!(stopwords, string.(['A':'z'; '0':'9'; '-'; '°']))
unique!(stopwords)

dir = "data/JAADtxt/"
files = dir .* readdir(dir)
@time txts = read.(files, String)
@time txts = [stem_all(Stemmer("english"), txt) for txt in txts]

crps = txts .|> StringDocument |> Corpus
remove_case!(crps)
@time update_lexicon!(crps)                  

[delete!(crps.lexicon, stopword) for stopword in stopwords];
@time for lexi in keys(crps.lexicon)
    if lexi ∉ gowords
        delete!(crps.lexicon, lexi)
    end
end

# _voca = crps.lexicon |> keys |> collect
# frq5 = _voca[(crps.lexicon |> values |> collect) .≤ 5]
# [delete!(crps.lexicon, trash) for trash in frq5];
voca = crps.lexicon |> keys |> collect |> sort

@time spX = crps |> DocumentTermMatrix |> dtm
@time _adjM = spX'spX
_adjM = sparse((_adjM - diagm(diag(_adjM)))) 
adjM = sparse(_adjM .> 100)

degree = sum(adjM, dims = 2) |> vec
histogram(degree, xlabel = "degree of node", ylabel = "frequency", size = [400, 400], xlims = [0, 100], label = :none, bins = 0:10:100, color = :black)
scatterhist(degree, xlabel = "degree of node", ylabel = "frequency", size = [400, 400], label = :none, color = :black, scale = :log10)
scatter(exp10.(0.25:0.25:4)[1:(end-1)], [count(exp10.(0:0.25:4)[k] .≤ degree .< exp10.(0:0.25:4)[k+1]) for k in 1:15], xscale = :log10, label = :none, size = [400, 400], title = "threshold = 10")

vocadegree = DataFrame(; voca, degree)
sort!(vocadegree, :degree, rev = true)
CSV.write("vocadegree.csv", vocadegree, bom = true)
AN = WordNetwork(adjM, voca)
sizeof(AN)
jldsave("G:/ElderlyAD/atopynetwork_cached.jld2"; AN)

# cp("G:/atopynetwork_cached.jld2", "atopynetwork_cached.jld2", force = true)


using Graphs
G = SimpleGraph(adjM)
unique(length.(connected_components(G)))