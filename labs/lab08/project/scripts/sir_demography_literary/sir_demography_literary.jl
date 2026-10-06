using DrWatson
@quickactivate "project"

include(srcdir("sir_model_demography.jl"))

using Random
using StatsPlots
using CSV
using Dates

tmax = 40.0

u0 = [
    990,
    10,
    0
]

p = [
    0.05,   # β - вероятность заражения
    10.0,   # c - частота контактов
    0.25,   # γ - интенсивность выздоровления
    0.01,   # μ - интенсивность смерти
    10.0    # λ_birth - интенсивность рождения
]

Random.seed!(1234)

demography_model = MakeSIRModel(
    u0,
    p
)

activate(
    demography_model
)

sir_run(
    demography_model,
    tmax
)

data_demography = out(
    demography_model
)

println(
    data_demography
)

@df data_demography plot(
    :t,
    [:S :I :R :D],
    label = ["S" "I" "R" "D"],
    xlab = "Время",
    ylab = "Численность",
    title = "SIR модель с демографическими событиями",
)

savefig(
    plotsdir(
        "sir_demography.png"
    )
)

mkpath(
    datadir("sims")
)

filename =
    "sir_demography_$(u0[1])_$(u0[2])_$(p[1])_$(p[2])_$(p[3])_$(p[4])_$(p[5]).csv"

CSV.write(
    datadir(
        "sims",
        filename
    ),
    data_demography
)

println(
    "Результаты сохранены в: ",
    datadir(
        "sims",
        filename
    )
)

println(
    "График сохранён в: ",
    plotsdir(
        "sir_demography.png"
    )
)
