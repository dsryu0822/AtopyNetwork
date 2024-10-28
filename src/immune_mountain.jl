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

wds = 1:5
folds = 1:5
pts = 90:1:100

N = length(txts)
# result = DataFrame(fd = Int64[], wd = Int64[], pt = Float64[], cover_score = Float64[])
# for dr = Base.product(folds, wds, pts, 0)
#     push!(result, dr)
# end; sort!(result, [:fd, :wd, :pt])

result = DataFrame(fd = Int64[], wd = Int64[], pt = Float64[], cp = Float64[], es = Float64[], cs1 = Float64[], cs2 = Float64[])
for wd = wds
    @load "/home/$(ENV["LOGNAME"])/temp/test_$wd.jld2"
    # W_

    wgtM_t_ = []
    θt = zeros(length(folds), length(pts))
    for fold = folds
        idx_v = fold:5:N; idx_t = setdiff(1:N, idx_v)
        wgtM_t = sum(W_[idx_t])
        push!(wgtM_t_, wgtM_t)
        nz_ = wgtM_t.nzval
        # nz_ = collect(wgtM_t_[fold][.!iszero.(wgtM_t_[fold])])
        θt[fold, :] .= percentile.(Ref(nz_), pts)
    end
    θv = zeros(N, length(pts))
    for iv = 1:N
        wgtM_v = W_[iv]
        nz_ = wgtM_v.nzval
        θv[iv, :] .= percentile.(Ref(nz_), pts)
    end

    for fold = folds
        idx_v = fold:5:N; idx_t = setdiff(1:N, idx_v)
        for ptk = eachindex(pts)
            adjM_t = wgtM_t_[fold] .> θt[fold, ptk]

            proper, dA, vol = zeros(N), zeros(N), zeros(N)
            @threads for iv = idx_v
                wgtM_v = W_[iv]
                adjM_v = wgtM_v .> θv[iv, ptk]
                proper[iv] = adjM_t != adjM_v)
                dA[iv] = count(adjM_t - adjM_v .< 0))
                vol[iv] = count(adjM_v)/2)
            end
            cp = 100count(iszero.(dA)) / length(dA)
            es = sum(vol) / length(vol)
            cs1 = cp * es
            cs2 = (vol'iszero.(dA)) / length(dA)

            try
                push!(result, [fold, wd, pts[ptk], cp, es, cs1, cs2])
                CSV.write("cover_score $device.csv", result)
            catch
                @error "error in $fold, $wd, $(pts[ptk])"
            end
        end
    end
end