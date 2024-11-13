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

wds = parse(Int64, last(device)):3:15
folds = 1:1
pts = 5:5:95

N = length(txts)
result = DataFrame(fd = Int64[], wd = Int64[], pt = Float64[], cp = Float64[], es = Float64[], cs1 = Float64[], cs2 = Float64[])
for wd = wds
    W_ = jldopen("../../temp/test_$wd.jld2")["W_"]
    # wd = 1; W_ = jldopen("weight_immune/test_$wd.jld2")["W_"]
    # @load "/home/$(ENV["LOGNAME"])/temp/test_$wd.jld2"

    wgtM_t_ = []
    θt = zeros(length(folds), length(pts))
    for fold = folds
        idx_v = fold:5:N; idx_t = setdiff(1:N, idx_v)
        wgtM_t = sum(W_[idx_t])
        push!(wgtM_t_, wgtM_t)
        nz_ = wgtM_t.nzval
        # θt[fold, :] .= percentile.(Ref(nz_), pts)
        θt[fold, :] .= pts' * maximum(nz_)
    end
    θv = zeros(N, length(pts))
    for iv = 1:N
        wgtM_v = W_[iv]
        nz_ = wgtM_v.nzval
        # θv[iv, :] .= percentile.(Ref(nz_), pts)
        θv[iv, :] .= pts' * maximum(nz_)
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
            cp = count(iszero.(dA)) / length(dA)
            es = sum(vol) / length(vol)
            cs1 = cp * es
            cs2 = (vol'iszero.(dA)) / length(dA)

            try
                push!(result, [fold, wd, pts[ptk], cp, es, cs1, cs2])
                CSV.write("cover_score percent $device [5, 95].csv", result)
            catch
                @error "error in $fold, $wd, $(pts[ptk])"
            end
        end
    end
end


# θ_ = sort(unique(wgtM_t_[1].nzval))
# n_node = []
# n_link = []
# @showprogress for θ = θ_
#     A = wgtM_t_[1] .> θ
#     push!(n_node, count(.!iszero.(sum(A, dims = 2))))
#     push!(n_link, sum(A))
# end
# pop!(n_node); pop!(n_link); pop!(θ_)

# plot(θ_, n_node, label = "n_node", lw = 2, scale = :log10)
# plot(θ_, n_link, label = "n_link", lw = 2, scale = :log10)

# heatmap(pts, 1:5, log10.(θt))
# heatmap(pts, 1:N, log10.(θv))