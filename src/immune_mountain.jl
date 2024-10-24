
include("../core/analyzer.jl")
cd("G:/corpus/"); pwd() # cd(@__DIR__)
using Base.Threads
@load "data_immune/cached_data.jld2"

folds = 1:1
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