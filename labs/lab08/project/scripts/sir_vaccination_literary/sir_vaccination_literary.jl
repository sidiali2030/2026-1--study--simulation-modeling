using DrWatson
@quickactivate "project"

include(srcdir("sir_model_vaccination.jl"))

using Random
using StatsPlots
using DataFrames

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

vaccination_time = 10.0

vaccination_fraction = 0.50

Random.seed!(1234)

model = MakeSIRModel(
    u0,
    p
)

activate(
    model
)

activate_vaccination!(
    model,
    vaccination_time,
    vaccination_fraction
)

sir_run(
    model,
    tmax
)

data_vaccination = out(
    model
)

println()

println(
    "========================================"
)

println(
    "РЕЗУЛЬТАТЫ"
)

println(
    "========================================"
)

println(
    "Время вакцинации = ",
    vaccination_time
)

println(
    "Доля вакцинации = ",
    vaccination_fraction * 100,
    "%"
)

println(
    "Максимальное число инфицированных = ",
    maximum(data_vaccination.I)
)

println(
    "Конечное S = ",
    data_vaccination.S[end]
)

println(
    "Конечное I = ",
    data_vaccination.I[end]
)

println(
    "Конечное R = ",
    data_vaccination.R[end]
)

plt = @df data_vaccination plot(
    :t,
    [:S :I :R],
    label = ["S" "I" "R"],
    xlabel = "Время",
    ylabel = "Численность",
    title = "SIR модель с вакцинацией",
    linewidth = 2,
    size = (900, 600)
)

vline!(
    plt,
    [vaccination_time],
    linestyle = :dash,
    label = "Вакцинация"
)

savefig(
    plt,
    plotsdir("sir_vaccination.png")
)

println()

println(
    "График сохранён: ",
    plotsdir("sir_vaccination.png")
)
