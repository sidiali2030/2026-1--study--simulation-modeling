using DrWatson
@quickactivate "project"

include(srcdir("sir_model.jl"))

using Random
using StatsPlots
using BenchmarkTools
using CSV
using Dates

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

Random.seed!(1234)

des_model = MakeSIRModel(
    u0,
    p
)

activate(
    des_model
)

sir_run(
    des_model,
    tmax
)

data_des = out(
    des_model
)

plt = @df data_des plot(
    :t,
    [:S :I :R],
    label = ["S" "I" "R"],
    xlab = "Время",
    ylab = "Численность",
    title = "Дискретно-событийная SIR модель",
)

savefig(
    plt,
    plotsdir("sir_des.png")
)

println(
    "График сохранён в: ",
    plotsdir("sir_des.png")
)

mkpath(
    datadir("sims")
)

filename =
    "sir_$(u0[1])_$(u0[2])_$(p[1])_$(p[2])_$(p[3]).csv"

CSV.write(
    datadir(
        "sims",
        filename
    ),
    data_des
)

println(
    "Результаты сохранены в: ",
    datadir(
        "sims",
        filename
    )
)
