using DrWatson
@quickactivate "project"

include(srcdir("sir_model_demography.jl"))

using Random
using StatsPlots
using CSV
using Dates

# Параметры модели
tmax = 40.0

u0 = [990, 10, 0] # S, I, R

# β, c, γ, μ, λ_birth
p = [
    0.05,   # β - вероятность заражения
    10.0,   # c - частота контактов
    0.25,   # γ - интенсивность выздоровления
    0.01,   # μ - интенсивность смерти
    10.0    # λ_birth - интенсивность рождения
]

Random.seed!(1234)

# Запуск модели
demography_model = MakeSIRModel(u0, p)

activate(demography_model)

sir_run(demography_model, tmax)

data_demography = out(demography_model)

# Вывод результатов
println(data_demography)

# Визуализация
@df data_demography plot(
    :t,
    [:S :I :R :D],
    labels = ["S" "I" "R" "D"],
    xlab = "Время",
    ylab = "Численность",
    title = "SIR модель с демографическими событиями",
)

savefig(
    plotsdir("sir_demography.png")
)

# Сохранение результатов в CSV
mkpath(datadir("sims"))

filename =
    "sir_demography_$(u0[1])_$(u0[2])_$(p[1])_$(p[2])_$(p[3])_$(p[4])_$(p[5]).csv"

CSV.write(
    datadir("sims", filename),
    data_demography
)

println(
    "Результаты сохранены в: ",
    datadir("sims", filename)
)

println(
    "График сохранён в: ",
    plotsdir("sir_demography.png")
)
