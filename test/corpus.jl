include("analyzer.jl")
packages = [:TextAnalysis, :Snowball, :CSV, :Graphs, :GraphRecipes, :LinearAlgebra, :NetworkLayout]
@showprogress for package in packages
    @eval using $(package)
end
mm = Plots.mm

rule_concat = ["skin barrier" => "skinxbarrier", "lamellar granul" => "lamellarxgranul", "antimicrobial peptid" => "antimicrobialxpeptid"]
rule_DNA = r"[ATCGU-]{4,}" => "xDNAseq"
rule_ref = r"\[\d*\]" => ' '
rule_numbers = string.(0:9) .=> ["O","l","Z","E","y","S","b","L","B","q"]
rule_ascii = setdiff(['!':'@'; '[':'`'; '{':'~'; '’'], ['-']) .=> ' '

rule_synonym = [" older " => " elder ", "antimicrobialxpeptid" => "AMP", "filaggrin" => "flg", "toll-likereceptors" => "TLRs"]

number2 = collect(Base.product(
    ["o","l","z","e","y","s","b","l","b","q"],
    ["o","l","z","e","y","s","b","l","b","q"]))
number2 = vec(prod.(number2))

# cd("//155.230.155.221/ty/Elderly AD/")
cd("G:/ElderlyAD/")

stopwords = CSV.read("./clinical-concepts/clinical-stopwords.txt", DataFrame)[:, 1]
append!(stopwords, string.(['A':'z'; '0':'9'; '-'; '°']))
append!(stopwords, number2)
unique!(stopwords)

dir = "data/JAADJID/"
files = dir .* readdir(dir)
# filter!(file -> !occursin.("melanoma", file), files)
@time txts = reduce.(*, readlines.(files))

txts = replace.(txts, rule_concat...)
txts = replace.(txts, rule_DNA)
txts = replace.(txts, rule_ref)
txts = replace.(txts, rule_numbers...)
txts = replace.(txts, rule_ascii...)

@time txts = [stem_all(Stemmer("english"), txt) for txt in txts]
txts = replace.(txts, rule_synonym...)

# count(occursin.("flg", txts))
# count(occursin.("age", txts) .&& occursin.("flg", txts))
# CSV.write("Flg appear.csv", DataFrame(title = files[occursin.("flg", txts)]))

crps = txts .|> StringDocument |> Corpus
remove_case!(crps)
@time update_lexicon!(crps)                  

[delete!(crps.lexicon, stopword) for stopword in stopwords];

_voca = crps.lexicon |> keys |> collect
frq1 = _voca[(crps.lexicon |> values |> collect) .≤ 5]
[delete!(crps.lexicon, trash) for trash in frq1];
voca = crps.lexicon |> keys |> collect |> sort
# CSV.write("voca.csv", DataFrame(; voca), bom = true)

# @time M_coo = CooMatrix(crps)
# @time sM_coo = coom(M_coo)
# M_coo |> propertynames
@time spX = crps |> DocumentTermMatrix |> dtm
@time adjM = spX'spX
# spX uses only 1.35% of the space

degree = sum(adjM, dims = 2) |> vec
# histogram(log10.(degree), legend = :none, xlabel = "log10(degree)", ylabel = "frequency", title = "Degree Distribution of Corpus Network")

vocadegree = DataFrame(; voca, degree)
sort!(vocadegree, :degree, rev = true)
CSV.write("vocadegree.csv", vocadegree, bom = true)
AN = WordNetwork(adjM, voca)
sizeof(AN)
jldsave("G:/ElderlyAD/atopynetwork_cached.jld2"; AN)

# cp("G:/atopynetwork_cached.jld2", "atopynetwork_cached.jld2", force = true)


# AN = []
# if isfile("G:/atopynetwork_cached.jld2")
#     jldopen("G:/atopynetwork_cached.jld2") do file
#         push!(AN, file["AN"])
#     end
# end
# AN = only(AN)
# vocadegree = CSV.read("//155.230.155.221/ty/Elderly AD/request5/vocadegree.csv", DataFrame)

selected = [
"age", "elder",
"skinxbarri", "lipid",
"pathogen", "environment", "allergen", "cytokin",
"keratinocyt", "corneocyt",
"kallikrein",
"flg", "amp", # "tlrs",
# "bailey", "volum", # "atlas", "beal",
"dermi", "syndrom", "peptid", "lifetim", "amino", "lineag"
]

node_id_ = findall(AN.nodes .∈ Ref(selected)); AN.nodes[node_id_]
toprate = 100 * [findfirst(vocadegree.voca .== select) for select in selected] ./ nrow(vocadegree)
deg = [vocadegree.degree[findfirst(vocadegree.voca .== select)] for select in selected]

percentile_ = []
# node_id = node_id_[1]
for node_id in node_id_
    neighbors = neighbor(AN, node_id)
    bit_stranger = .!(selected .∈ Ref(neighbors.node))
    if any(bit_stranger)
        neighbors = DataFrame([Matrix(neighbors); [node_id_[bit_stranger] selected[bit_stranger] fill(0, count(bit_stranger))]], names(neighbors))
    end
    rank_selected = findall(neighbors.node .∈ Ref(selected))

    neighbors_selected = neighbors[rank_selected, :]
    percentile = rank_selected ./ nrow(neighbors)
    neighbors_selected[!, "percentil"] = percentile
    sort!(neighbors_selected, :id)
    push!(percentile_, percentile)
end
pctM = stack(percentile_)
pctM = (pctM + pctM') / 2
adjM = (pctM .< 0.10)
adjM[diagind(adjM)] .= 0



_selected = replace(selected,
"flg" => "Flg", "tlrs" => "TLRs", "amp" => "AMP", "skinxbarri" => "skin barrier",
"keratinocyt" => "keratinocyte", "corneocyt" => "corneocyte", "cytokin" => "cytokine",
"peptid" => "peptide",
reverse.(rule_concat)...)

using Random
Random.seed!(27)
x_, y_ = eachrow(0.5 .+ 4rand(2, length(_selected)))
# x_, y_ = eachrow(stackstress(adjM))
x_[1:2] .= [4, 4.3]
y_[1:2] .= [4, 4.5]

x_[3:4] .= [0.9, 0.7]
y_[3:4] .= [4.3, 3.8]
x_[9:10] .= [1.7, 1.23]
y_[9:10] .= [4.1, 3.2]
x_[12] = 0.8
y_[12] = 3.6


x_[5:7] .= [3.1, 3.2, 2.9]
y_[5:7] .= [1.1, 2.2, 1.7]

# x_[1:14] += 0.1randn(14)
# y_[1:14] += 0.1randn(14)
x_[15:end] += 4cospi.((15:length(_selected)) ./ 3)
y_[15:end] += 4sinpi.((15:length(_selected)) ./ 3)

# graphplot(adjM[1:14, 1:14], x = x_[1:14], y = y_[1:14], names = _selected[1:14], la = 0.5, mc = :white, ms = 0)
nodearg = (; ms = 20, msw = 0, alpha = 0.5, shape = :rect)
q0 = graphplot(adjM, x = x_, y = y_, la = 0.1, ms = 0)
q1 = scatter(q0, x_, y_, text = Plots.text.(_selected, :gray), ms = 0)
plot!(q1, lims = [0, 6], size = [600,600])

q2 = graphplot!(q0, adjM[1:14, 1:14], x = x_[1:14], y = y_[1:14], la = 0.5, mc = :white, ms = 0)
scatter!(q2, x_[15:end], y_[15:end], text = Plots.text.(_selected[15:end], :gray), ms = 0)
scatter!(q2, x_[1:2], y_[1:2], text = Plots.text.(_selected[1:2], :black), color = 1; nodearg...)
scatter!(q2, x_[[3,4,9,10,12]], y_[[3,4,9,10,12]], text = Plots.text.(_selected[[3,4,9,10,12]], :black), color = 2; nodearg...)
scatter!(q2, x_[5:7], y_[5:7], text = Plots.text.(_selected[5:7], :black), color = 3; nodearg...)
scatter!(q2, x_[[8,11,13,14]], y_[[8,11,13,14]], text = Plots.text.(_selected[[8,11,13,14]], :black); color = :gray, nodearg...)
plot!(q2, lims = [0, 6], size = [600,600])

plot(q1, framestyle = :box, dpi = 300); png("q1")
plot(q2, framestyle = :box, dpi = 300); png("q2")

plot(q1, q2, layout = (1,2), size = [1200, 600], margin = 5mm, framestyle = :box, dpi = 300)

N = nrow(vocadegree)
npcrank = rand(1:N, 1000)
npc = vocadegree[npcrank,:]

colors_ = [1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3]
shapes_ = repeat([:dtriangle, :utriangle, :diamond, :rect, :pentagon, :hexagon, :circle], 2)

p4_1 = plot(legend = :bottomright, size = [660, 660],
xscale = :log10, yscale = :log10, xlims = 10.0 .^ [4,8], ylims = 10.0 .^ [-1, 2],
yformatter = y -> "$(y) %", xlabel = "Degree of node", ylabel = "Top percentage of node",
xticks = 10.0 .^ [4,6,8], yticks = 10.0 .^ [-1, 0, 1, 2], yflip = true)
scatter!(p4_1, npc.degree, 100*(npcrank ./ N), color = :black, label = :none, alpha = 0.1, shape = :xcross)
[scatter!(p4_1, deg[[k]], toprate[[k]], label = _selected[k],
mc = :white, ma = 0.5, ms = 10, msw = 2, msc = colors_[k], shape = shapes_[k]) for k in 1:14]
p4_1

p4_2 = plot(legend = :bottomright, size = [660, 660],
ylims = [-5, 100],
yformatter = y -> "$(y) %", xlabel = "Degree of node", ylabel = "Top percentage of node",
yticks = [0,50,100], yflip = true)
scatter!(p4_2, npc.degree, 100*(npcrank ./ N), color = :black, label = :none, alpha = 0.5, shape = :xcross)
[scatter!(p4_2, deg[[k]], toprate[[k]], label = _selected[k],
mc = :white, ma = 0.5, ms = 10, msw = 2, msc = colors_[k], shape = shapes_[k]) for k in 1:14]
p4_2

plot(
    plot(p4_2, xlabel = "", legend = :none)
  , plot(p4_1, legend = :topleft)
  , framestyle = :box
  , layout = (2,1), size = [500, 900], leftmargin = 7mm, rightmargin = 5mm, grid = false, dpi = 300)
png("deg_toppct")

BC = betweenness_centrality(Graph(adjM))
CC = closeness_centrality(Graph(adjM))
EC = eigenvector_centrality(Graph(adjM))

w_age = []
w_elder = []
for u in selected
    _neighbor = neighbor(AN, u)
    push!(w_age, _neighbor.weight[findfirst(_neighbor.node .== "age")])
    push!(w_elder, _neighbor.weight[findfirst(_neighbor.node .== "elder")])
end

CSV.write("query_table.csv", DataFrame(; selected, deg, toprate, w_age, w_elder))

ANsymmetric = deepcopy(AN.adj)
ANsymmetric[diagind(ANsymmetric)] .= 0
ANtemp = ANsymmetric[1:100:end, 1:100:end]
colors_ = [[1,1,2,2,2,2,2,3,3,3,:gray,:gray,:gray,:gray]; fill(:white, 327)]
alphas_ = 0.5*(colors_ .!= :white)
@time graphplot(ANtemp,
la = 0.01, method = :spectral,
shape = :rect, ma = alphas_, msw = 0, msa = 0,
nodecolor = colors_)
plot!(size = [600,600], dpi = 300)
png("overview")



println("")
for node ∈ selected
    CSV.write("$(node).csv", neighbor(AN, "$node"), bom = true)
end
neighbor(AN, "ad")
neighbor(AN, "keratinocyt")
neighbor(AN, "lipid")

# test = DataFrame(spX, voca)[:, [:age, :flg, :amp, :keratinocyt, :corneocyt, :cytokin, :peptid, :dermi, :syndrom, :lipid, :pathogen, :environment, :allergen, :kallikrein, :lifetim, :amino, :lineag]]
