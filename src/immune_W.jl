include("../core/analyzer.jl")
if Sys.iswindows()
    cd("G:/corpus/")
else
    cd("/home/$(ENV["LOGNAME"])/g/corpus/")
end; @info pwd()
using Base.Threads, Dates
device = gethostname()
@info "$(now()) - $device $(nthreads()) threads"

@load "data_immune/cached_data.jld2"
N = length(txts)
crps_ = [Corpus(crps[[k]]) for k in 1:N]
for _crps in crps_
    _crps.lexicon = crps.lexicon
end

for wd = 1:15
    @info "Window size: $wd"
    W_ = Array{SparseMatrixCSC}(undef, N)
    @showprogress @threads for c in 1:N
        W_[c] = CooMatrix(crps_[c], window = wd, normalize = false).coom
    end
    @save "cached_$(wd).jld2" W_
end