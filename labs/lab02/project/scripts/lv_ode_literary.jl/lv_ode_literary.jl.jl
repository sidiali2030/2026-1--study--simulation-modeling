using DrWatson
@quickactivate "project"

using DifferentialEquations
using DataFrames
using StatsPlots
using LaTeXStrings
using Plots
using Statistics
using FFTW

script_name = splitext(basename(PROGRAM_FILE))[1]

mkpath(plotsdir(script_name))
mkpath(datadir(script_name))

function lotka_volterra!(du, u, p, t)
    x, y = u
    α, β, δ, γ = p

    @inbounds begin
        du[1] = α * x - β * x * y
        du[2] = δ * x * y - γ * y
    end

    nothing
end

p_lv = [
    0.1,  # α: скорость размножения жертв
    0.02, # β: скорость поедания жертв хищниками
    0.01, # δ: коэффициент конверсии пищи в хищников
    0.3,  # γ: смертность хищников
]

u0_lv = [40.0, 9.0]

tspan_lv = (0.0, 200.0)
dt_lv = 0.01

prob_lv = ODEProblem(
    lotka_volterra!,
    u0_lv,
    tspan_lv,
    p_lv,
)

sol_lv = solve(
    prob_lv,
    Tsit5();
    dt = dt_lv,
    reltol = 1e-8,
    abstol = 1e-10,
    saveat = 0.1,
    dense = true,
)

df_lv = DataFrame()

df_lv[!, :t] = sol_lv.t
df_lv[!, :prey] = [u[1] for u in sol_lv.u]
df_lv[!, :predator] = [u[2] for u in sol_lv.u]

df_lv[!, :dprey_dt] =
    p_lv[1] .* df_lv.prey .-
    p_lv[2] .* df_lv.prey .* df_lv.predator

df_lv[!, :dpredator_dt] =
    p_lv[3] .* df_lv.prey .* df_lv.predator .-
    p_lv[4] .* df_lv.predator

println("="^60)
println("Модель Лотки-Вольтерры (хищник-жертва)")
println("="^60)

println("\nПараметры модели:")
println("α (скорость размножения жертв) = ", p_lv[1])
println("β (скорость поедания жертв) = ", p_lv[2])
println("δ (коэффициент конверсии) = ", p_lv[3])
println("γ (смертность хищников) = ", p_lv[4])

println("\nНачальные условия:")
println("Жертвы (x0) = ", u0_lv[1])
println("Хищники (y0) = ", u0_lv[2])

x_star = p_lv[4] / p_lv[3]
y_star = p_lv[1] / p_lv[2]

println("\nСтационарные точки (положения равновесия):")
println("x* = γ/δ = ", round(x_star, digits = 3))
println("y* = α/β = ", round(y_star, digits = 3))

plt1 = plot(
    df_lv.t,
    [df_lv.prey df_lv.predator],
    label = [L"Жертвы (x)" L"Хищники (y)"],
    xlabel = "Время",
    ylabel = "Популяция",
    title = "Модель Лотки-Вольтерры: Динамика популяций",
    linewidth = 2,
    legend = :topright,
    grid = true,
    size = (900, 500),
    color = [:green :red],
)

hline!(
    plt1,
    [x_star],
    color = :green,
    linestyle = :dash,
    alpha = 0.5,
    label = "x* (равновесие жертв)",
)

hline!(
    plt1,
    [y_star],
    color = :red,
    linestyle = :dash,
    alpha = 0.5,
    label = "y* (равновесие хищников)",
)

plt2 = plot(
    df_lv.prey,
    df_lv.predator,
    label = "Фазовая траектория",
    xlabel = "Популяция жертв (x)",
    ylabel = "Популяция хищников (y)",
    title = "Фазовый портрет системы",
    color = :blue,
    linewidth = 1.5,
    grid = true,
    size = (800, 600),
    legend = :topright,
)

step = 50

for i in 1:step:length(df_lv.prey)-step
    plot!(
        plt2,
        [df_lv.prey[i], df_lv.prey[i+step]],
        [df_lv.predator[i], df_lv.predator[i+step]],
        arrow = :closed,
        color = :blue,
        alpha = 0.3,
        label = false,
    )
end

scatter!(
    plt2,
    [x_star],
    [y_star],
    color = :black,
    markersize = 8,
    label = "Стационарная точка (x*, y*)",
)

x_range =
    LinRange(0, maximum(df_lv.prey) * 1.1, 100)

y_nullcline =
    p_lv[1] ./ (p_lv[2] .* x_range)

plot!(
    plt2,
    x_range,
    y_nullcline,
    color = :red,
    linestyle = :dash,
    linewidth = 1.5,
    label = "Изоклина хищников (dy/dt=0)",
)

y_range =
    LinRange(0, maximum(df_lv.predator) * 1.1, 100)

x_nullcline =
    p_lv[4] ./
    (p_lv[3] .* ones(length(y_range)))

plot!(
    plt2,
    x_nullcline,
    y_range,
    color = :green,
    linestyle = :dash,
    linewidth = 1.5,
    label = "Изоклина жертв (dx/dt=0)",
)

plt3 = plot(
    df_lv.t,
    [df_lv.dprey_dt df_lv.dpredator_dt],
    label = [L"dx/dt" L"dy/dt"],
    xlabel = "Время",
    ylabel = "Скорость изменения",
    title = "Производные популяций",
    linewidth = 1.5,
    legend = :topright,
    grid = true,
    size = (900, 400),
    color = [:green :red],
)

hline!(
    plt3,
    [0],
    color = :black,
    linestyle = :solid,
    alpha = 0.3,
    label = false,
)

df_lv[!, :prey_pct_change] =
    df_lv.dprey_dt ./ df_lv.prey .* 100

df_lv[!, :predator_pct_change] =
    df_lv.dpredator_dt ./ df_lv.predator .* 100

plt4 = plot(
    df_lv.t,
    [
        df_lv.prey_pct_change
        df_lv.predator_pct_change
    ],
    label = [
        L"dx/dt / x (\%)"
        L"dy/dt / y (\%)"
    ],
    xlabel = "Время",
    ylabel = "Относительное изменение, %",
    title = "Относительные темпы роста",
    linewidth = 1.5,
    legend = :topright,
    grid = true,
    size = (900, 400),
    color = [:green :red],
)

function compute_fft(signal, dt)
    n = length(signal)

    spectrum = abs.(rfft(signal))

    freq = rfftfreq(n, 1 / dt)

    return freq, spectrum
end

freq_prey, spectrum_prey =
    compute_fft(
        df_lv.prey .- mean(df_lv.prey),
        dt_lv,
    )

freq_predator, spectrum_predator =
    compute_fft(
        df_lv.predator .- mean(df_lv.predator),
        dt_lv,
    )

plt5 = plot(
    freq_prey,
    [spectrum_prey spectrum_predator],
    label = [L"Жертвы (x)" L"Хищники (y)"],
    xlabel = "Частота",
    ylabel = "Амплитуда",
    title = "Спектральный анализ (Фурье)",
    linewidth = 1.5,
    xscale = :log10,
    yscale = :log10,
    legend = :topright,
    grid = true,
    size = (800, 400),
    color = [:green :red],
)

if length(spectrum_prey) > 0
    idx_prey =
        argmax(spectrum_prey[2:end]) + 1

    dominant_freq_prey =
        freq_prey[idx_prey]

    period_prey =
        1 / dominant_freq_prey

    println(
        "\nДоминирующая частота колебаний жертв: ",
        round(dominant_freq_prey, digits = 4),
        " Гц",
    )

    println(
        "Период колебаний жертв: ",
        round(period_prey, digits = 2),
        " единиц времени",
    )
end

plt6 = plot(
    layout = (3, 2),
    size = (1200, 900),
)

plot!(
    plt6[1],
    df_lv.t,
    df_lv.prey,
    label = L"x(t)",
    color = :green,
    linewidth = 2,
    title = "Популяция жертв",
    grid = true,
)

plot!(
    plt6[2],
    df_lv.t,
    df_lv.predator,
    label = L"y(t)",
    color = :red,
    linewidth = 2,
    title = "Популяция хищников",
    grid = true,
)

plot!(
    plt6[3],
    df_lv.prey,
    df_lv.predator,
    label = false,
    color = :blue,
    linewidth = 1.5,
    title = "Фазовый портрет",
    xlabel = L"x",
    ylabel = L"y",
    grid = true,
)

scatter!(
    plt6[3],
    [x_star],
    [y_star],
    color = :black,
    markersize = 5,
    label = "(x*, y*)",
)

plot!(
    plt6[4],
    df_lv.t,
    [df_lv.dprey_dt df_lv.dpredator_dt],
    label = [L"dx/dt" L"dy/dt"],
    color = [:green :red],
    linewidth = 1.5,
    title = "Скорости изменения",
    grid = true,
    legend = :topright,
)

plot!(
    plt6[5],
    freq_prey,
    spectrum_prey,
    label = L"x",
    color = :green,
    linewidth = 1.5,
    title = "Спектр жертв",
    xscale = :log10,
    yscale = :log10,
    grid = true,
)

plot!(
    plt6[6],
    df_lv.t,
    [
        df_lv.prey_pct_change
        df_lv.predator_pct_change
    ],
    label = [L"dx/x" L"dy/y"],
    color = [:green :red],
    linewidth = 1.5,
    title = "Относительные изменения",
    grid = true,
    legend = :topright,
)

println("\n" * "="^60)
println("Анализ результатов")
println("="^60)

println("\nОсновные статистики:")

println(
    "Жертвы: min = ",
    round(minimum(df_lv.prey), digits = 2),
    ", max = ",
    round(maximum(df_lv.prey), digits = 2),
    ", mean = ",
    round(mean(df_lv.prey), digits = 2),
)

println(
    "Хищники: min = ",
    round(minimum(df_lv.predator), digits = 2),
    ", max = ",
    round(maximum(df_lv.predator), digits = 2),
    ", mean = ",
    round(mean(df_lv.predator), digits = 2),
)
