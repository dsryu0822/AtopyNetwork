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
exception = Dict("T cell" => "T-cell", "B cell" => "B-cell")

# @time raw = read.(readdir("data_immune/test", join = true), String)
@time raw = [read.(readdir("data_immune/NI_article", join = true), String);
             read.(readdir("data_immune/NRI_review", join = true), String);
             read.(readdir("data_immune/NI_review", join = true), String)];
txts = replace.(raw, cleaning...)
txts = replace.(txts, exception...)
txts = replace.(txts, "T-cells" => "T-cell ")
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

# @save "D:/cached_data.jld2" txts crps voca

include("../core/analyzer.jl")
cd("G:/corpus/"); pwd() # cd(@__DIR__)
using Base.Threads
@load "data_immune/cached_data.jld2"

folds = 2:5
wds = 1:9
pts = 5:5:95

N = length(txts)
result = DataFrame(fd = [], wd = [], pt = [], cover_score = [])
for wd = wds
    crps_ = [Corpus(crps[[k]]) for k in 1:N]
    wgtM_ = []
    @time for _crps in crps_
        _crps.lexicon = crps.lexicon
        push!(wgtM_, CooMatrix(_crps, window = wd).coom)
    end
    for fold = folds
        idx_v = fold:5:N; idx_t = setdiff(1:N, idx_v)
        wgtM_t = sum(wgtM_[idx_t])
        nz_ = collect(wgtM_t[.!iszero.(wgtM_t)])
        @threads for pt = pts
            θ = percentile(nz_, pt)
            adjM_t = wgtM_t .> θ
            
            dA = []
            vol = []
            @showprogress for iv = idx_v
                wgtM_v = wgtM_[iv]
                adjM_v = wgtM_v .> θ
                push!(dA, count(adjM_t - adjM_v .< 0))
                push!(vol, count(adjM_v)/2)
            end
            cover_score = (vol'iszero.(dA)) / length(dA)
            # cover_score = 100count(iszero.(dA)) / length(dA)            

            try
                push!(result, [fold, wd, pt, cover_score])
                CSV.write("cover_score $fold.csv", result)
            catch
                @error "error in $fold, $wd, $pt"
            end
        end
    end
end

result = CSV.read("cover_score.csv", DataFrame)
result = result[result.fd .== 1, :]
sort!(result, [:wd, :pt])

sfc = reshape(result.cover_score, 11, :)
surface(1:9, 0:10:100, sfc, alpha = .5, zlims = [0, 100], xlabel = "window size", ylabel = "percentile", zlabel = "cover score")
png("surface_cover_score")

degree = vec(sum(adjM_t, dims = 2))
IN = WordNetwork(adjM_t, voca)
# @save String(@__DIR__) * "/../immune_netwrork.jld2" IN

@load "G:/corpus/immune_network/IN.jld2"
neighbor(IN, "cell")
voca = replace.(IN.nodes, reverse.(('0':'9') .=> ('０':'９'))...)
dg = vec(sum(IN.adj, dims = 2))
vd = DataFrame(; voca, dg); sort!(vd, :dg, rev = true)
logscatter(dg[.!iszero.(dg)], scale = :log10, size = [400, 400], color = :black, xlabel = L"d", ylabel = L"p(d)", label = :none)

for l = 1:5
    idx_sbgrph = findall(.!iszero.(dg .> 25))[l:5:end]
    G = SimpleGraph(IN.adj[idx_sbgrph, idx_sbgrph])
    gargs = (; size = 2*[2160, 2160], alpha = .5, nodesize = 0, la = .1, lims = [-.3, .3], fontsize = 21)
    @time p1 = graphplot(G, names = voca[idx_sbgrph]; gargs...)
    png(p1, "G:/corpus/immune_network/zoom_$l")
end

graphplot(G, names = voca[idx_sbgrph], lims = [-.3, .3])