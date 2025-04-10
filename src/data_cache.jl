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
allow_ascii = [48:57; 65:90; 97:122; 43; 45] # Char.(allow_ascii)

target = "data_celldeath"
files = vcat([joinpath.(root, files) for (root, dirs, files) = walkdir(target)]...)
raw = read.(files, String);
txts = replace.(raw, cleaning...)
txts = replace.(txts, cell...)
txts = replace.(txts, "cells" => "cell ")
@time txts = [stem_all(Stemmer("english"), txt) for txt in txts]
txts = replace.(txts, r" (\d)" => s"\1")
txts = replace.(txts, r"(\d) \+" => s"\1+")
txts = replace.(txts, r"[ATCGU-]{4,}" => ":xDNAseq")
txts = replace.(txts, r"([a-z])+(\d)+" => "")
txts = replace.(txts, (('0':'9') .=> ('０':'９'))...)

crps = txts .|> StringDocument |> Corpus
remove_case!(crps)
@time update_lexicon!(crps)
[delete!(crps.lexicon, stopword) for stopword in stopwords];
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
write("takealook.txt", join(voca, ", "))

# @save "cached_data.jld2" txts crps voca

N = length(txts)
crps_ = [Corpus(crps[[k]]) for k in 1:N]
for _crps in crps_
    _crps.lexicon = crps.lexicon
end

# try mkdir(joinpath("cached", target)) catch e end
# for wd = 1:10
#     @info "Window size: $wd"
#     W_ = Array{SparseMatrixCSC}(undef, N)
#     @showprogress @threads for c in 1:N
#         W_[c] = CooMatrix(crps_[c], window = wd, normalize = false).coom
#     end
#     @save joinpath("cached", target, "wd_$(wd).jld2") W_
# end

############################################
#                   test                   #
############################################
W_ = Array{SparseMatrixCSC}(undef, N)
@showprogress @threads for c in 1:N
    W_[c] = Int64.(CooMatrix(crps_[c], window = 1, normalize = false).coom)
end

W = sum(W_)
idcs = findall(!iszero, W)
df_edge = DataFrame(Source = [], Target = [], Weight = [])
for idx = idcs
    i, j = collect(idx.I)
    if i > j
        push!(df_edge, [voca[[i, j]]; W[idx]])
    end
end
sort!(df_edge, :w, rev = true)
CSV.write("$(target)_weight.csv", df_edge)


for year = string.(2010:2024)
    W = sum(W_[occursin.(year, files)])
    idcs = findall(!iszero, W)
    df_edge = DataFrame(Source = [], Target = [], Weight = [])
    for idx = idcs
        i, j = collect(idx.I)
        if i > j
            push!(df_edge, [voca[[i, j]]; W[idx]])
        end
    end
    sort!(df_edge, :w, rev = true)
    CSV.write("$(year)_weight.csv", df_edge)
end