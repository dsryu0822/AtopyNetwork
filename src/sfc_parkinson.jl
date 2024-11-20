include("../core/analyzer.jl")
if Sys.iswindows()
    cd("G:/corpus/")
else
    cd("/home/$(ENV["LOGNAME"])/g/corpus/")
end; @info pwd()
using Base.Threads, Dates
device = gethostname()
@info "$(now()) - $device $(nthreads()) threads"
cd("data_nature")

hargs = (; xlabel = "Window size", ylabel = "Percentile")

result = CSV.read("data_nature/cover_score [95, 100].csv", DataFrame)
wd_ = unique(result.wd)
pt_ = unique(result.pt)
gdf1 = groupby(result, :fd)
surface(wd_, 90:99, sum([reshape(df.cp, length(pt_), :) for df in gdf1]) ./ 5, zformatter = x -> "$(trunc(Int64, 100x)) %", yscale = :log10; hargs...)
