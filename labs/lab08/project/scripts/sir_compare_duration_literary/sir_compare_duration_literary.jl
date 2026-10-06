using DrWatson
@quickactivate "project"

using Random
using StatsPlots

include(srcdir("sir_model.jl"))
include(srcdir("sir_model_fixed.jl"))

tmax = 40.0

u0 = [
    990,
    10,
    0
]

p = [
    0.05,
    10.0,
    0.25
]

seed = 1234

Random.seed!(seed)

stoch_model = MakeSIRModel(
    u0,
    p
)

activate(
    stoch_model
)

sir_run(
    stoch_model,
    tmax
)

data_stoch = out(
    stoch_model
)

Random.seed!(seed)

fixed_model = MakeSIRModelFixed(
    u0,
    p
)

activate_fixed(
    fixed_model
)

sir_run_fixed(
    fixed_model,
    tmax
)

data_fixed = out_fixed(
    fixed_model
)

plt = plot(
    data_stoch.t,
    data_stoch.I,
    label = "I — случайная длительность",
    xlabel = "Время",
    ylabel = "Число инфицированных",
    title = "Сравнение длительности болезни",
    linewidth = 2,
)

plot!(
    plt,
    data_fixed.t,
    data_fixed.I,
    label = "I — фиксированная длительность",
    linewidth = 2,
    linestyle = :dash,
)

savefig(
    plt,
    plotsdir("sir_duration_comparison.png")
)

println("График сохранён:")
println(
    plotsdir("sir_duration_comparison.png")
)
