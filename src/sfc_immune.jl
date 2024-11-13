include("../core/analyzer.jl")
if Sys.iswindows()
    cd("G:/corpus/")
else
    cd("/home/$(ENV["LOGNAME"])/g/corpus/")
end; @info pwd()
using Base.Threads, Dates
device = gethostname()
@info "$(now()) - $device $(nthreads()) threads"

hargs = (; xlabel = "Window size", ylabel = "Percentile", size = [800, 800])

result = CSV.read("data_immune/cover_score percentile [90, 100].csv", DataFrame)
gdf1 = groupby(result, :fd)
surface(1:15, 90:1:100, sum([reshape(df.cp, 11, :) for df in gdf1]) ./ 5, zformatter = x -> "$(trunc(Int64, 100x)) %"; hargs...)
# surface(1:15, 90:1:100, sum([reshape(df.es, 11, :) for df in gdf1]) ./ 5; hargs...)
surface(1:15, 90:1:100, sum([reshape(df.cs1, 11, :) for df in gdf1]) ./ 5, camera = [30, 60]; hargs...)
surface(1:15, 90:1:100, sum([reshape(df.cs2, 11, :) for df in gdf1]) ./ 5; hargs...)

scatter(
    vcat([df.wd for df in groupby(result, :wd)]...),
    vcat([df.cs2 for df in groupby(result, :wd)]...),
    ms = 2, msw = 0, label = :none, color = :black, xticks = 1:15, xlabel = "Window size"
)
scatter(
    vcat([df.pt for df in groupby(result, :pt)]...),
    vcat([df.cs2 for df in groupby(result, :pt)]...),
    ms = 2, msw = 0, label = :none, color = :black, xticks = 90:100, xlabel = "Percentile"
)
png("temp")

gdf1_2 = groupby(result[result.wd .== 1, :], :fd)
p1 = plot()
for df in gdf1_2 plot!(p1, df.cs2) end

result = CSV.read("cover_score percentile [5, 95].csv", DataFrame)
gdf2 = groupby(result, :fd)
heatmap(1:15, 5:5:95, sum([reshape(df.cp, 19, :) for df in gdf2]) ./ 5; hargs...)
# surface(1:15, 5:5:95, sum([reshape(df.es, 19, :) for df in gdf2]) ./ 5; hargs...)
heatmap(1:15, 5:5:95, sum([reshape(df.cs1, 19, :) for df in gdf2]) ./ 5; hargs...)
heatmap(1:15, 5:5:95, sum([reshape(df.cs2, 19, :) for df in gdf2]) ./ 5; hargs...)

gdf3 = groupby(result[result.pt .> 40, :], :fd)
surface(1:15, 45:5:95, sum([reshape(df.cs2, 11, :) for df in gdf3]) ./ 5; hargs...)
