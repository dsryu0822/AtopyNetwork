include("analyzer.jl")
packages = [:TextAnalysis, :Snowball, :CSV]
for package in ProgressBar(packages)
    @eval using $(package)
end

rule_concat = "skin disease" => "skinxdisease"
rule_DNA = r"[ATCGU-]{4,}" => "xDNAseq"
rule_ref = r"\[\d*\]" => ' '
rule_numbers = string.(0:9) .=> ["O","l","Z","E","y","S","b","L","B","q"]
rule_ascii = setdiff(['!':'@'; '[':'`'; '{':'~'; '’'], ['-']) .=> ' '

rule_synonym = " older " => " elder "

number2 = collect(Base.product(
    ["o","l","z","e","y","s","b","l","b","q"],
    ["o","l","z","e","y","s","b","l","b","q"]))
number2 = vec(prod.(number2))

cd("//155.230.155.221/ty/Elderly AD/")
# cd("G:/Elderly AD/")

stopwords = CSV.read("./clinical-concepts/clinical-stopwords.txt", DataFrame)[:, 1]
append!(stopwords, string.(['A':'z'; '0':'9'; '-'; '°']))
append!(stopwords, number2)
unique!(stopwords)

dir = "data/JAAD/"
files = dir .* readdir(dir)
# filter!(file -> !occursin.("melanoma", file), files)
@time txts = reduce.(*, readlines.(files))

txts = replace.(txts, rule_concat)
txts = replace.(txts, rule_DNA)
txts = replace.(txts, rule_ref)
txts = replace.(txts, rule_numbers...)
txts = replace.(txts, rule_ascii...)

@time txts = [stem_all(Stemmer("english"), txt) for txt in txts]
txts = replace.(txts, rule_synonym)

crps = txts .|> StringDocument |> Corpus
remove_case!(crps)
@time update_lexicon!(crps)                  

[delete!(crps.lexicon, stopword) for stopword in stopwords];

_voca = crps.lexicon |> keys |> collect
frq1 = _voca[(crps.lexicon |> values |> collect) .≤ 5]
[delete!(crps.lexicon, trash) for trash in frq1];
voca = crps.lexicon |> keys |> collect |> sort
# CSV.write("voca.csv", DataFrame(; voca), bom = true)

@time spX = crps |> DocumentTermMatrix |> dtm
@time adjM = spX'spX                         
# spX uses only 1.35% of the space

degree = sum(adjM, dims = 2) |> vec
# histogram(log10.(degree), legend = :none, xlabel = "log10(degree)", ylabel = "frequency", title = "Degree Distribution of Corpus Network")

vocadegree = DataFrame(; voca, degree)
sort!(vocadegree, :degree, rev = true)
CSV.write("vocadegree.csv", vocadegree, bom = true)

atopynetwork = WordNetwork(adjM, voca)
jldsave("atopynetwork_cached.jld2"; atopynetwork)
cp("atopynetwork_cached.jld2", "G:/atopynetwork_cached.jld2")

# for node in ["age", "psoriasi", "skin", "eczema", "ad", "elder", "older", "skinxdiseas"]
#     plot_subgraph(AN, node, 5, 2); png(node)
#     CSV.write("$node.csv", neighbor(AN, node), bom = true)
# end
# neighbor(AN, "skinxdiseas")
# plot_subgraph(AN, "skinxdiseas", 5, 2); png("skinxdiseas")
# neighbor(AN, "skin")
# CSV.write("AD.csv", neighbor(AN, "AD"))
# CSV.write("eczema.csv", neighbor(AN, "eczema"))
# CSV.write("elder.csv", neighbor(AN, "elder"))
# CSV.write("older.csv", neighbor(AN, "older"))

["age", "psoriasi", "skin", "eczema", "ad", "elder", "skinxdiseas"]