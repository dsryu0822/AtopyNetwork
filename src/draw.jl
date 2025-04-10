using Plots, CSV, DataFrames

sbjt = "immune"
sbjt = "parkinson"
sbjt = "celldeath"
sbjt = "genetics"
sbjt = "alzheimer"
sbjt = "rheumatoid"
cd("G:/corpus/data_$sbjt")
file = CSV.read("coverability.csv", DataFrame)

wd_ = unique(file.wd)
pt_ = unique(file.pt)
cp__ = sum(reshape(file.cp, length(pt_), 5, length(wd_)), dims = 2)[:,1,:] / 5

hm = heatmap(wd_, eachindex(pt_), cp__, color = :grays)
sf = surface(wd_, eachindex(pt_), cp__, clims = (0, 1), zlims = [0, 1])
plot(hm, sf, size = [1200, 500], plot_title = "$sbjt")
png("cp")

