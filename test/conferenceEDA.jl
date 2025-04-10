using Base.Threads: @threads, nthreads # Base.Threads.nthreads()
using ProgressMeter
packages = [:CSV, :DataFrames, :Plots, :DecisionTree, :GLM]
@showprogress @threads for package in packages
    @eval using $(package)
end

cd("G:/network")
default(alpha = .5, lw = 0)
# data = vcat(CSV.read.(["$year.csv" for year in 2022:2024], DataFrame)...)
data = CSV.read("graph.csv", DataFrame)
feature = names(data)[3:end]

data0 = data[.!data.y, :]
data1 = data[data.y, :]

hargs = (; bins = 100, normalize = :pdf)
for ftr = feature
    p = plot(title = ftr)
    histogram!(p, data0[:, ftr], label = "reject"; hargs...)
    histogram!(p, data1[:, ftr], label = "accept"; hargs...)
    png(p, "EDA/$ftr.png")
end

to_take_log = ["diameter", "max_core", "max_weight", "min_cut", "num_leaf", "radius"]
for ftr = to_take_log
    p = plot(title = "log_$ftr")
    histogram!(p, log1p.(data0[:, ftr]), label = "reject"; hargs...)
    histogram!(p, log1p.(data1[:, ftr]), label = "accept"; hargs...)
    png(p, "EDA/$ftr.png")
end
data[:, to_take_log] .= log10.(data[:, to_take_log])
default()

count(.!data.y) / nrow(data)

out = glm(@formula(y ~ year + pareto + max_central + max_weight + max_core + max_ind + max_pagerank + min_cut + min_dom + diameter + radius + num_leaf + coef_clustering + coef_assortativity + entropy + mad), data, Bernoulli(), LogitLink())
prob = GLM.predict(out, data[:, [:y, :year, :pareto, :max_central, :max_weight, :max_core, :max_ind, :max_pagerank, :min_cut, :min_dom, :diameter, :radius, :num_leaf, :coef_clustering, :coef_assortativity, :entropy, :mad]])
θ = 0.5
count(data.y .== (prob .> θ)) / nrow(data)

acc_ = []
for θ = 0:0.01:1
    push!(acc_, count(data.y .== (prob .> θ)) / nrow(data))
end
plot(0:0.01:1, acc_)

model = DecisionTreeClassifier(max_depth=10)
features = Matrix(data[:, [:year, :pareto, :max_central, :max_weight, :max_core, :max_ind, :max_pagerank, :min_cut, :min_dom, :diameter, :radius, :num_leaf, :coef_clustering, :coef_assortativity, :entropy, :mad]]);
labels = data.y
DecisionTree.fit!(model, features, labels)
print_tree(model)
sum(labels .== DecisionTree.predict(model, features)) / length(labels)
