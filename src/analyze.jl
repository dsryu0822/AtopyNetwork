using Base.Threads, Dates
device = gethostname()
sbjt = first(ARGS)
@info "$(now()) - $device $(nthreads()) threads - $sbjt"

include("../core/analyzer.jl")
if Sys.iswindows()
    cd("G:/corpus/data_$sbjt")
else
    cd("/home/$(ENV["LOGNAME"])/g/corpus/data_$sbjt")
end; @info pwd()

# wds = parse(Int64, last(device)):2:10
wds = 1:10
folds = 1:5
pts = sort(100*(1 .- logrange(1e-4, .05, 10)))

result = DataFrame(fd = Int64[], wd = Int64[], pt = Float64[], cp = Float64[], es = Float64[], cs1 = Float64[], cs2 = Float64[])
for wd = wds
    if Sys.iswindows()
        W_ = jldopen("../_cached/$sbjt/cached_$wd.jld2")["W_"]
    else
        W_ = jldopen("/home/$(ENV["LOGNAME"])/_cached/$sbjt/cached_$wd.jld2")["W_"]
    end
    N = length(W_)

    wgtM_t_ = []
    θt = zeros(length(folds), length(pts))
    for fold = folds
        idx_v = fold:5:N; idx_t = setdiff(1:N, idx_v)
        wgtM_t = sum(W_[idx_t])
        push!(wgtM_t_, wgtM_t)
        nz_ = wgtM_t.nzval
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
            adjM_t = wgtM_t_[fold] .≥ θt[fold, ptk]

            proper, dA, vol = zeros(N), zeros(N), zeros(N)
            @threads for iv = idx_v
                wgtM_v = W_[iv]
                adjM_v = wgtM_v .≥ θv[iv, ptk]
                proper[iv] = adjM_t != adjM_v
                dA[iv] = count(adjM_t - adjM_v .< 0)
                vol[iv] = count(adjM_v)/2
            end
            proper, dA, vol = proper[idx_v], dA[idx_v], vol[idx_v]
            cp = count(iszero.(dA)) / length(dA)
            es = sum(vol) / length(vol)
            cs1 = cp * es
            cs2 = (vol'iszero.(dA)) / length(dA)

            try
                push!(result, [fold, wd, pts[ptk], cp, es, cs1, cs2])
                CSV.write("$(device)_coverability.csv", result)
            catch
                @error "error in $fold, $wd, $(pts[ptk])"
            end
        end
    end
end