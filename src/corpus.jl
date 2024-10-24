include("../core/analyzer.jl")

cd("G:/ElderlyAD/data")
txts_ = String[]
N_ = []
for dir = readdir()
    files = readdir(dir, join = true)
    @time txts = [stem_all(Stemmer("english"), txt) for txt in read.(files, String)]
    push!(N_, length(txts))
    append!(txts_, txts)
end; cd("G:/ElderlyAD")
journal = UnitRange.([1; (cumsum(N_) .+ 1)[1:(end-1)]], cumsum(N_))

crps = txts_ .|> StringDocument |> Corpus
remove_case!(crps)
@time update_lexicon!(crps)

stopwords = [stem(Stemmer("english"), word) for word in lowercase.(CSV.read("./clinical-concepts/clinical-stopwords.txt", DataFrame)[:, 1])]
unique!(stopwords)
[delete!(crps.lexicon, stopword) for stopword in stopwords];

gowords = [stem(Stemmer("english"), word) for word in lowercase.(CSV.read("gowords.csv", DataFrame)[:, 1])]
unique!(gowords)
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
CSV.write("dataed/df_sem.csv", DataFrame(spX, voca), bom = true)

####################################################
#                                                  #
#                     vanilla                      #
#                                                  #
####################################################
@time _adjM = spX'spX
_adjM = sparse((_adjM - diagm(diag(_adjM))))
adjM = sparse(_adjM .> 0)
# adjM = _adjM

degree = sum(adjM, dims = 2) |> vec
histogram(degree, xlabel = "degree of node", ylabel = "frequency", size = [400, 400], label = :none, color = :black)
scatterhist(degree, xlabel = "degree of node", ylabel = "frequency", size = [400, 400], label = :none, color = :black, scale = :log10)
scatter(exp10.(0.25:0.25:4)[1:(end-1)], [count(exp10.(0:0.25:4)[k] .≤ degree .< exp10.(0:0.25:4)[k+1]) for k in 1:15], xscale = :log10, label = :none, size = [400, 400], title = "threshold = 10")

vocadegree = DataFrame(; voca, degree)
sort!(vocadegree, :degree, rev = true)
CSV.write("dataed/vanilla/vocadegree.csv", vocadegree, bom = true)
AN = WordNetwork(adjM, voca)
sizeof(AN)

####################################################
#                                                  #
#                   multilayer                     #
#                                                  #
####################################################
layer = []
for jn in journal
    _adjM = spX[jn, :]'spX[jn, :]
    _adjM = sparse((_adjM - diagm(diag(_adjM))))
    push!(layer, sparse(_adjM .> 0))
end
adjM = sum(layer)
degree = sum(adjM, dims = 2) |> vec

vocadegree = DataFrame(; voca, degree)
sort!(vocadegree, :degree, rev = true)
CSV.write("dataed/multilayer/vocadegree.csv", vocadegree, bom = true)
AN = WordNetwork(adjM, voca)
CSV.write("dataed/multilayer/N(age).csv", neighbor(AN, "age"), bom = true)

####################################################
#                                                  #
#                     pressed                      #
#                                                  #
####################################################
adjM = sparse(sum(layer) .== 4)
degree = sum(adjM, dims = 2) |> vec
vocadegree = DataFrame(; voca, degree)
sort!(vocadegree, :degree, rev = true)
CSV.write("dataed/pressed/vocadegree.csv", vocadegree, bom = true)
AN = WordNetwork(adjM, voca)
CSV.write("dataed/pressed/N(atopi).csv", neighbor(AN, "atopi"), bom = true)


scatterhist(degree[degree .> 0], xlabel = "degree of node", ylabel = "frequency", size = [400, 400],
label = :none, color = :black, yscale = :log10, ticks = exp10.(0:1:3))

####################################################
#                                                  #
#      Co occurrence matrix one by one corpus      #
#                                                  #
####################################################
@time C = CooMatrix(crps, normalize = false)
n_isolated_ = []
for θ = 0:15
    push!(n_isolated_, sum(iszero.(sum(coom(C) .> θ, dims = 2))))
end
scatter(0:15, n_isolated_, title = "number of isolated nodes", xlabel = "threshold", xticks = 0:15)
adjM = coom(C) .> 5

degree = sum(adjM, dims = 2) |> vec; sum(degree)/2length(degree)
vocadegree = DataFrame(; voca, degree)
sort!(vocadegree, :degree, rev = true)
CSV.write("dataed/coom/vocadegree.csv", vocadegree, bom = true)
AN = WordNetwork(adjM, voca)
CSV.write("dataed/coom/N(atopi).csv", neighbor(AN, "atopi"), bom = true)
CSV.write("dataed/coom/N(age).csv", neighbor(AN, "age"), bom = true)

neighbor(AN, "apposit")

histogram(degree, xlabel = "degree of node", ylabel = "frequency", size = [400, 400], label = :none, color = :black)
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
logscatter(degree, scale = :log10, color = :black, size = [400, 400])

using Graphs
# G = SimpleGraph(adjM[.!iszero.(degree), .!iszero.(degree)])
G = SimpleGraph(adjM)
id = Dict(voca .=> eachindex(voca))

path_ = []
push!(path_, yen_k_shortest_paths(G, id["atopi"], id["age"]).paths[1])
voca[path_[1]]
push!(path_, yen_k_shortest_paths(G, id["keratinocyt"], id["age"]).paths[1])
voca[path_[2]]
push!(path_, yen_k_shortest_paths(G, id["keratinocyt"], id["atopi"]).paths[1])
voca[path_[3]]
