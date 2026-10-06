using DrWatson
@quickactivate "project"

include(srcdir("seir_model.jl"))

using Random
using StatsPlots
using DataFrames

tmax = 40.0

u0 = [
    990,
    0,
    10,
    0
]

p = [
    0.05,
    10.0,
    0.5,
    0.25
]

println("========================================")
println("SEIR МОДЕЛЬ")
println("========================================")

println("Начальное состояние:")
println("S0 = ", u0[1])
println("E0 = ", u0[2])
println("I0 = ", u0[3])
println("R0 = ", u0[4])

println()

println("Параметры:")
println("β = ", p[1])
println("c = ", p[2])
println("σ = ", p[3])
println("γ = ", p[4])

println(
    "Средний латентный период = ",
    1 / p[3]
)

println(
    "Средняя длительность болезни = ",
    1 / p[4]
)

Random.seed!(1234)

seir_model = MakeSEIRModel(
    u0,
    p
)

activate(seir_model)

seir_run(
    seir_model,
    tmax
)

data_seir = out(seir_model)

peak_I = maximum(data_seir.I)

peak_index = argmax(data_seir.I)

peak_time = data_seir.t[peak_index]

final_S = data_seir.S[end]
final_E = data_seir.E[end]
final_I = data_seir.I[end]
final_R = data_seir.R[end]

println()

println("========================================")
println("РЕЗУЛЬТАТЫ")
println("========================================")

println(
    "Максимальное число инфицированных = ",
    peak_I
)

println(
    "Время пика = ",
    round(
        peak_time,
        digits = 3
    )
)

println(
    "Конечное S = ",
    final_S
)

println(
    "Конечное E = ",
    final_E
)

println(
    "Конечное I = ",
    final_I
)

println(
    "Конечное R = ",
    final_R
)

println(
    "Проверка населения = ",
    final_S +
    final_E +
    final_I +
    final_R
)

plt = @df data_seir plot(
    :t,
    [:S :E :I :R],
    label = ["S" "E" "I" "R"],
    xlabel = "Время",
    ylabel = "Численность",
    title = "Дискретно-событийная SEIR модель",
    linewidth = 2,
    size = (900, 600)
)

savefig(
    plt,
    plotsdir("seir_des.png")
)

println()

println(
    "График сохранён: ",
    plotsdir("seir_des.png")
)
