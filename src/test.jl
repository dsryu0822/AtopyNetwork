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
cell = Dict("T cell" => "T-cell", "B cell" => "B-cell")

cd("data_test")
@time raw = [read.(filter(x -> occursin(".txt", x), readdir("NI_article_fig", join = true)), String);
             read.(filter(x -> occursin(".txt", x), readdir("NI_review_fig", join = true)), String);
             read.(filter(x -> occursin(".txt", x), readdir("NRI_fig", join = true)), String)];
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


# default(fontfamily = "")
reduced_term = collect(keys(crps.lexicon))
exportdir = [readdir("NI_article_fig", join = true); readdir("NI_review_fig", join = true); readdir("NRI_fig", join = true)]
filter!(x -> occursin(".txt", x), exportdir)

th = 40
mkdir.(["NI_article_png", "NI_article_csv", "NI_review_png", "NI_review_csv", "NRI_png", "NRI_csv"])
csv_ = replace.(exportdir, "fig" => "csv", ".txt" => ".csv")
png_ = replace.(exportdir, "fig" => "png", ".txt" => ".png")
@showprogress for ck = 1:29
    try
        CM = CooMatrix(crps[ck], reduced_term, window = 1, normalize = false)
        W = CM.coom
        A = W .≥ th
        idx_node = .!iszero.(vec(sum(A, dims = 2)))
        rA = A[idx_node, idx_node]
        terms = CM.terms[idx_node]
        terms = replace.(terms, (('０':'９') .=> ('0':'9'))...)
        G = SimpleGraph(rA)

        temp = [terms DataFrame(Int64.(rA), terms)]
        rename!(temp, :x1 => :terms)
        CSV.write(csv_[ck], temp)

        graphplot(G, names = terms, size = [1600, 1600], ms = 0, la = .1)
        png(png_[ck])
    catch
        @error "error in $ck"
    end
end