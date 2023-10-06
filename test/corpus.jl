using ProgressBars
packages = [:TextAnalysis, :DataFrames, :Plots, :CSV, :Snowball]
for package in ProgressBar(packages)
    @eval using $(package)
end
stopwords = CSV.read("stopwords.csv", DataFrame)
stopwords = stopwords.voca[ismissing.(stopwords.except)]
append!(stopwords, string.(['A':'z'; '0':'9'; '-'; '°']))

# rule_numbers = string.(0:9) .=> ["０","１","２","３","４","５","６","７","８","９"]
rule_DNA = r"[ATCGU-]{4,}" => "xDNAseq"
rule_ref = r"\[\d*\]" => ' '
rule_numbers = string.(0:9) .=> ["O","l","Z","E","y","S","b","L","B","q"]
rule_ascii = setdiff(['!':'@'; '[':'`'; '{':'~'; '’'], ['-']) .=> ' '

files = "data/" .* readdir("data")
txts = reduce.(*, readlines.(files))

txts = replace.(txts, rule_DNA)
txts = replace.(txts, rule_ref)
txts = replace.(txts, rule_numbers...)
txts = replace.(txts, rule_ascii...)

@time txts = [stem_all(Stemmer("english"), txt) for txt in txts] # 8.740254 seconds

crps = txts .|> StringDocument |> Corpus
@time update_lexicon!(crps)                   # 8.698965 seconds (41.99 M allocations: 2.263 GiB, 4.27% gc time)

[delete!(crps.lexicon, stopword) for stopword in stopwords];
voca = crps |> lexicon |> keys |> collect |> sort
# CSV.write("voca.csv", DataFrame(; voca), bom = true)

@time spX = crps |> DocumentTermMatrix |> dtm # 8.528082 seconds (42.46 M allocations: 2.295 GiB, 3.32% gc time, 0.09% compilation time)
@time adjM = spX'spX                          # 9.169203 seconds (19 allocations: 12.511 GiB, 0.07% gc time)
# spX uses only 1.35% of the space

degree = sum(adjM, dims = 2) |> vec
histogram(log10.(degree), legend = :none, xlabel = "log10(degree)", ylabel = "frequency", title = "Degree Distribution of Corpus Network")

vocafreq = DataFrame(; voca, degree)
sort!(vocafreq, :degree, rev = true)
CSV.write("vocafreq.csv", vocafreq, bom = true)


idx_atop = findfirst(voca .== "atop")
atop = DataFrame(voc = voca[adjM[idx_atop, :] .> 0], weight = adjM[idx_atop, adjM[idx_atop, :] .> 0])
sort!(atop, :weight, rev = true)

idx_IL_13 = findfirst(voca .== "IL-lE")
IL_13 = DataFrame(voc = voca[adjM[idx_IL_13, :] .> 0], weight = adjM[idx_IL_13, adjM[idx_IL_13, :] .> 0])
sort!(IL_13, :weight, rev = true)
