include("../core/analyzer.jl")
cd("G:/corpus/"); pwd() # cd(@__DIR__)
using Base.Threads

# stopwords = [stem(Stemmer("english"), word) for word in lowercase.(CSV.read("./clinical-concepts/clinical-stopwords.txt", DataFrame)[:, 1])]
# unique!(stopwords)
# gowords = [stem(Stemmer("english"), word) for word in lowercase.(CSV.read("gowords.csv", DataFrame)[:, 1])]
# unique!(gowords)
# @save "G:/corpus/words.jld2" stopwords gowords
@load "words.jld2" 
gowords = [gowords; titlecase.(gowords)]
stopwords = [stopwords; titlecase.(stopwords)]

cleaning = Dict(" " => " ", "None\r\n" => "", "-" => " ", "+" => " ")
cell = Dict("T cell" => "T-cell", "B cell" => "B-cell", "NK cell" => "NK-cell")

cd("data_nature");
# @time raw = read.(readdir("Parkinson", join = true), String);
@time raw = read.(readdir("Cell_Death_Disease", join = true), String);
# cd("data_immune"); @time raw = [
#     read.(readdir("NI_article", join = true), String);
#     read.(readdir("NRI_review", join = true), String);
#     read.(readdir("NI_review", join = true), String)];
txts = replace.(raw, cleaning...)
txts = replace.(txts, cell...)
txts = replace.(txts, "cells" => "cell ")
@time txts = [stem_all(Stemmer("english"), txt) for txt in txts]
txts = replace.(txts, r" (\d)" => s"\1")
txts = replace.(txts, r"(\d) \+" => s"\1+")
txts = replace.(txts, r"[ATCGU-]{4,}" => ":xDNAseq")
txts = replace.(txts, r"([a-z])+(\d)+" => "")
txts = replace.(txts, (('0':'9') .=> ('０':'９'))...)
# @save "data_immune/raw.jld2" raw
# @load "data_immune/txts_cached.jld2"

crps = txts .|> StringDocument |> Corpus
remove_case!(crps) # 370763 -> 370330
@time update_lexicon!(crps) # 4576 doc in 104s
[delete!(crps.lexicon, stopword) for stopword in stopwords];
allow_ascii = [48:57; 65:90; 97:122; 43; 45; 65296:65305]
@time for (k, v) in crps.lexicon
    _lexi = collect(k)
    try
        if !all(codepoint.(_lexi) .∈ Ref(allow_ascii)) ||
        !(islowercase(first(k)) || isuppercase(first(k))) ||
        v ≤ 1

            delete!(crps.lexicon, k)
        end
    catch
        @error "error in $v, $k"
        delete!(crps.lexicon, k)
    end
end
voca = replace.(crps.lexicon |> keys |> collect, reverse.(('0':'9') .=> ('０':'９'))...)
open("takealook.txt", "w") do io
    println(io, voca)
end

@save "cached_data.jld2" txts crps voca

N = length(txts)
crps_ = [Corpus(crps[[k]]) for k in 1:N]
for _crps in crps_
    _crps.lexicon = crps.lexicon
end

for wd = 1:10
    @info "Window size: $wd"
    W_ = Array{SparseMatrixCSC}(undef, N)
    @showprogress @threads for c in 1:N
        W_[c] = CooMatrix(crps_[c], window = wd, normalize = false).coom
    end
    @save "cached_wd_$(wd).jld2" W_
end
# degree = vec(sum(adjM_t, dims = 2))
# IN = WordNetwork(adjM_t, voca)
# # @save String(@__DIR__) * "/../immune_netwrork.jld2" IN

# @load "G:/corpus/immune_network/IN.jld2"
# neighbor(IN, "cell")
# voca = replace.(IN.nodes, reverse.(('0':'9') .=> ('０':'９'))...)
# dg = vec(sum(IN.adj, dims = 2))
# vd = DataFrame(; voca, dg); sort!(vd, :dg, rev = true)
# logscatter(dg[.!iszero.(dg)], scale = :log10, size = [400, 400], color = :black, xlabel = L"d", ylabel = L"p(d)", label = :none)

# for l = 1:5
#     idx_sbgrph = findall(.!iszero.(dg .> 25))[l:5:end]
#     G = SimpleGraph(IN.adj[idx_sbgrph, idx_sbgrph])
#     gargs = (; size = 2*[2160, 2160], alpha = .5, nodesize = 0, la = .1, lims = [-.3, .3], fontsize = 21)
#     @time p1 = graphplot(G, names = voca[idx_sbgrph]; gargs...)
#     png(p1, "G:/corpus/immune_network/zoom_$l")
# end

# graphplot(G, names = voca[idx_sbgrph], lims = [-.3, .3])

# roulette = rand(1:20, length(txts))

# n_term = zeros(Int64, 20)
# @showprogress @threads for rl = 1:20
#     crps = txts[roulette .≤ rl] .|> StringDocument |> Corpus
#     remove_case!(crps) # 370763 -> 370330
#     update_lexicon!(crps) # 4576 doc in 104s
#     [delete!(crps.lexicon, stopword) for stopword in stopwords];
#     allow_ascii = [48:57; 65:90; 97:122; 43; 45; 65296:65305]
#     for (k, v) in crps.lexicon
#         _lexi = collect(k)
#         try
#             if !all(codepoint.(_lexi) .∈ Ref(allow_ascii)) ||
#             !(islowercase(first(k)) || isuppercase(first(k))) ||
#             v ≤ 1

#                 delete!(crps.lexicon, k)
#             end
#         catch
#             @error "error in $v, $k"
#             delete!(crps.lexicon, k)
#         end
#     end
#     n_term[rl] = length(crps.lexicon)
# end

# n_paper = [count(roulette .≤ rl) for rl in 1:20]
# plot(n_paper, n_term, legend = :none, size = [500, 500], color = :black, lw = 2,
# xticks = [410, 8826], yticks = [18634, 134934], yformatter = :plain)