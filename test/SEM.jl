using StructuralEquationModels

observed_vars = [:x1, :x2, :x3, :y1, :y2, :y3, :y4, :y5, :y6, :y7, :y8]
latent_vars = [:ind60, :dem60, :dem65]

graph = @StenoGraph begin

    # loadings
    ind60 → fixed(1)*x1 + x2 + x3
    dem60 → fixed(1)*y1 + y2 + y3 + y4
    dem65 → fixed(1)*y5 + y6 + y7 + y8

    # latent regressions
    ind60 → dem60
    dem60 → dem65
    ind60 → dem65

    # variances
    _(observed_vars) ↔ _(observed_vars)
    _(latent_vars) ↔ _(latent_vars)

    # covariances
    y1 ↔ y5
    y2 ↔ y4 + y6
    y3 ↔ y7
    y8 ↔ y4 + y6

end

partable = ParameterTable(
    latent_vars = latent_vars,
    observed_vars = observed_vars,
    graph = graph)

data = example_data("political_democracy")
model = Sem(
    specification = partable,
    data = data
)
model_fit = sem_fit(model)
fit_measures(model_fit)
sem_summary(model_fit)
update_estimate!(partable, model_fit)
update_partable!(partable, model_fit, se_hessian(model_fit), :se)
sem_summary(partable)

using DataFrames, Distributions

result = DataFrame(partable)[:, [:from, :parameter_type, :to, :estimate, :se]]
result[!, :z] = abs.(result.estimate ./ result.se)
result[!, :p] = 2(1 .- cdf.(Normal(), abs.(result.z)))

result

using CSV
test = CSV.read("test.csv", DataFrame)

observed_vars = Symbol.(names(test))
# latent_vars = [:C1, :C2, :C3, :C4]
latent_vars = [:C1]

graph = @StenoGraph begin
    # loadings
    # C1 → fixed(1)*age + flg + lifetim + amino + lineag
    # C2 → fixed(1)*amp + keratinocyt + syndrom + lipid + environment
    # C3 → fixed(1)*corneocyt + cytokin + peptid
    # C4 → fixed(1)*dermi + pathogen + allergen + kallikrein

    # # latent regressions
    # C1 → C2
    # C2 → C4
    # C3 → C4

    # variances
    _(observed_vars) ↔ _(observed_vars)
    # _(latent_vars) ↔ _(latent_vars)

    # covariances
    age → flg + lipid
    amino → lineag
    # lineag → amino
    keratinocyt → syndrom
    corneocyt → cytokin
    peptid ↔ dermi + pathogen + allergen + kallikrein
end

partable = ParameterTable(
    latent_vars = latent_vars,
    observed_vars = observed_vars,
    graph = graph)
model = Sem(
    specification = partable,
    data = test
)
model_fit = sem_fit(model)
fit_measures(model_fit)
sem_summary(model_fit)
update_estimate!(partable, model_fit)
update_partable!(partable, model_fit, se_hessian(model_fit), :se)
sem_summary(partable)

result = DataFrame(partable)[:, [:from, :parameter_type, :to, :estimate, :se]]
result[!, :z] = abs.(result.estimate ./ result.se)
result[!, :p] = 2(1 .- cdf.(Normal(), abs.(result.z)))

result



observed_vars = [:age, :dermi, :cytokin]
latent_vars = [:C1, :C2]
graph = @StenoGraph begin
    C1 → age + cytokin
    C2 → dermi
    
    # variances
    _(observed_vars) ↔ _(observed_vars)
    _(latent_vars) ↔ _(latent_vars)
    
    # covariances
    age ↔ dermi
end

partable = ParameterTable(
    latent_vars = latent_vars,
    observed_vars = observed_vars,
    graph = graph)
model = Sem(
    specification = partable,
    data = test
)
