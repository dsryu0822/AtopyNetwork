celldeath_ = first.(CSV.read.("셀데스/" .* string.(2010:2024) .* "_weight.csv", DataFrame), 500)

_nodes = sort(vcat(celldeath_...), :Weight, rev = true)
nodes = unique(reshape(vcat(reshape.([_nodes.Source, _nodes.Target], 1, :)...), :))
idcs = 1:length(nodes)
dni = Dict(nodes .=> idcs)

_x = (idcs .+ 10) .* cos.(idcs)
_y = (idcs .+ 10) .* sin.(idcs)

p1 = plot(size = [2000, 2000], lims = [-100, 100], legend = :none, framestyle = :none)

for year = 2010:2024
    clldth = celldeath_[year - 2009]
    u_ = [dni[ky] for ky in clldth.Source]
    v_ = [dni[ky] for ky in clldth.Target]
    clldth.u_ = u_
    clldth.v_ = v_
    clldth.w_ = 10log10.(clldth.Weight) ./ maximum(log10.(clldth.Weight))

    p2 = deepcopy(p1)
    for dr = eachrow(clldth)
        plot!(p2, [_x[dr.u_], _x[dr.v_]], [_y[dr.u_], _y[dr.v_]], c = :black, lw = dr.w_, alpha = .1)
    end
    p3 = scatter(p2, _x, _y, txt = nodes, msw = 0, mc = :white, ma = .5, ms = 20, shape = :rect)
    png(p3, "$(year)_weight.png")
end

for year = 2010:2024
    clldth = celldeath_[year - 2009]
    clldth.Weight .= log10.(clldth.Weight)
    rename!(clldth, :Weight => "Weight_$year")
end
panel = outerjoin(celldeath_..., on = [:Source, :Target])

Plots.plotly()
p4 = plot(size = [1900, 1200])
for i in 1:100
    plot!(p4, collect(panel[i, 3:end]) ./ panel[i, 3], label = "$(panel.Source[i]) - $(panel.Target[i])")
end
p4
savefig(p4, "링크 트렌드")