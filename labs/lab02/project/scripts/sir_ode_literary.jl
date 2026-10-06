# # Модель SIR
#
# **Цель:** исследовать динамику распространения инфекционного заболевания
# с помощью трёхпараметрической модели SIR.
#
# В модели население разделяется на три группы:
#
# - `S(t)` — восприимчивые;
# - `I(t)` — инфицированные;
# - `R(t)` — выздоровевшие.
#
# В трёхпараметрической форме используются параметры:
#
# - `β` — вероятность передачи инфекции при контакте;
# - `c` — среднее число контактов;
# - `γ` — скорость выздоровления.
#
# Система дифференциальных уравнений имеет вид:
#
# $$
# \frac{dS}{dt} = -\beta c \frac{IS}{N},
# $$
#
# $$
# \frac{dI}{dt} = \beta c \frac{IS}{N} - \gamma I,
# $$
#
# $$
# \frac{dR}{dt} = \gamma I,
# $$
#
# где
#
# $$
# N = S + I + R.
# $$

# ## Инициализация проекта и загрузка пакетов
#
# Активируем окружение DrWatson и подключаем библиотеки,
# необходимые для решения системы ОДУ, обработки данных,
# визуализации и оценки производительности.

using DrWatson
@quickactivate "project"

using DifferentialEquations
using SimpleDiffEq
using Tables
using DataFrames
using StatsPlots
using LaTeXStrings
using Plots
using BenchmarkTools

# Имя текущего скрипта используется для создания отдельных каталогов
# с графиками и результатами вычислений.

script_name = splitext(basename(PROGRAM_FILE))[1]

mkpath(plotsdir(script_name))
mkpath(datadir(script_name))

# ## Определение модели
#
# Функция `sir_ode!` задаёт правую часть системы дифференциальных
# уравнений SIR.
#
# Вектор состояния имеет вид:
#
# $$
# u = (S, I, R),
# $$
#
# а вектор параметров:
#
# $$
# p = (\beta, c, \gamma).
# $$

function sir_ode!(du, u, p, t)
    (S, I, R) = u
    (β, c, γ) = p

    N = S + I + R

    @inbounds begin
        du[1] = -β * c * I / N * S
        du[2] = β * c * I / N * S - γ * I
        du[3] = γ * I
    end

    nothing
end

# ## Параметры модели
#
# Зададим шаг интегрирования, продолжительность моделирования,
# начальное состояние популяции и параметры распространения инфекции.

δt = 0.1
tmax = 40.0
tspan = (0.0, tmax)

# Начальные условия:
#
# - 990 восприимчивых;
# - 10 инфицированных;
# - 0 выздоровевших.

u0 = [990.0, 10.0, 0.0]

# Параметры модели:
#
# - β = 0.05;
# - c = 10;
# - γ = 0.25.

p = [0.05, 10.0, 0.25]

# ## Базовое репродуктивное число
#
# Для трёхпараметрической модели SIR базовое репродуктивное число
# определяется выражением:
#
# $$
# R_0 = \frac{c\beta}{\gamma}.
# $$
#
# Если `R₀ > 1`, число инфицированных на начальном этапе может расти.
# Если `R₀ < 1`, эпидемический процесс затухает.

R0 = (p[2] * p[1]) / p[3]

# ## Решение системы дифференциальных уравнений
#
# Создаём объект `ODEProblem`, после чего численно решаем систему
# на заданном временном интервале.

prob_ode = ODEProblem(sir_ode!, u0, tspan, p)

sol_ode = solve(prob_ode, dt = δt)

# ## Подготовка результатов
#
# Результат решения преобразуем в `DataFrame`.
# Добавим время и общую численность популяции.

df_ode = DataFrame(Tables.table(sol_ode'))

rename!(df_ode, ["S", "I", "R"])

df_ode[!, :t] = sol_ode.t

df_ode[!, :N] =
    df_ode.S + df_ode.I + df_ode.R

# ## Вывод параметров модели
#
# Выведем используемые параметры, значение базового репродуктивного
# числа и начальные условия.

println("Параметры модели SIR:")
println("β (вероятность заражения) = ", p[1])
println("c (среднее число контактов) = ", p[2])
println("γ (скорость выздоровления) = ", p[3])
println("R0 = c * β / γ = ", round(R0, digits = 3))

println(
    "Средняя продолжительность болезни = ",
    round(1 / p[3], digits = 2),
    " дней",
)

println(
    "Начальные условия: S0 = ",
    u0[1],
    ", I0 = ",
    u0[2],
    ", R0 = ",
    u0[3],
)

# ## Динамика основных групп населения
#
# Построим основной график модели, содержащий одновременно
# траектории `S(t)`, `I(t)` и `R(t)`.

plt1 = @df df_ode plot(
    :t,
    [:S :I :R],
    label = [L"S(t)" L"I(t)" L"R(t)"],
    xlabel = "Время, дни",
    ylabel = "Количество людей",
    title = "Модель SIR: Динамика эпидемии",
    linewidth = 2,
    legend = :right,
    grid = true,
    size = (800, 500),
)

# На график добавим основные параметры модели.

annotate!(
    plt1,
    maximum(df_ode.t) * 0.7,
    maximum(df_ode.N) * 0.8,
    text(
        "Параметры:\nβ = $(p[1])\nc = $(p[2])\nγ = $(p[3])\nR0 = $(round(R0, digits=2))",
        8,
        :left,
    ),
)

# ## Динамика инфицированных
#
# Отдельно рассмотрим изменение числа инфицированных `I(t)`.

plt2 = @df df_ode plot(
    :t,
    :I,
    label = L"I(t)",
    xlabel = "Время, дни",
    ylabel = "Количество инфицированных",
    title = "Динамика числа зараженных",
    color = :red,
    linewidth = 2,
    fill = (0, 0.3, :red),
    grid = true,
    size = (800, 400),
)

# ### Пик эпидемии
#
# Найдём максимальное число инфицированных и момент времени,
# когда достигается этот максимум.

peak_idx = argmax(df_ode.I)
peak_time = df_ode.t[peak_idx]
peak_value = df_ode.I[peak_idx]

vline!(
    plt2,
    [peak_time],
    color = :black,
    linestyle = :dash,
    label = false,
    linewidth = 1,
)

annotate!(
    plt2,
    peak_time,
    peak_value * 1.05,
    text(
        "Пик: $(round(peak_value, digits=1)) на $(round(peak_time, digits=1)) день",
        8,
        :top,
    ),
)

# ## Логарифмический масштаб
#
# Логарифмическая шкала позволяет лучше исследовать начальный
# экспоненциальный рост числа инфицированных.

plt3 = @df df_ode plot(
    :t,
    :I,
    label = L"I(t)",
    xlabel = "Время, дни",
    ylabel = "Количество инфицированных (лог. масштаб)",
    title = "Экспоненциальный рост (лог. шкала)",
    yscale = :log10,
    color = :red,
    linewidth = 2,
    grid = true,
    size = (800, 400),
)

# ## Доли групп населения
#
# Представим состояния `S`, `I` и `R` не в абсолютных значениях,
# а в процентах от общей численности населения.

plt4 = @df df_ode plot(
    :t,
    [:S :I :R] ./ df_ode.N .* 100,
    label = [L"S(t)/N" L"I(t)/N" L"R(t)/N"],
    xlabel = "Время, дни",
    ylabel = "Доля популяции, %",
    title = "Динамика эпидемии (в процентах)",
    linewidth = 2,
    legend = :right,
    grid = true,
    size = (800, 500),
)

# ### Порог коллективного иммунитета
#
# Для `R₀ > 1` теоретический порог коллективного иммунитета:
#
# $$
# H = 1-\frac{1}{R_0}.
# $$

if R0 > 1
    herd_immunity_threshold = (1 - 1 / R0) * 100

    hline!(
        plt4,
        [herd_immunity_threshold],
        color = :purple,
        linestyle = :dash,
        label =
            "Порог коллективного иммунитета ($(round(herd_immunity_threshold, digits=1))%)",
        linewidth = 1.5,
    )
end

# ## Фазовый портрет
#
# Фазовый портрет показывает взаимосвязь между числом
# восприимчивых `S` и числом инфицированных `I`.

plt5 = plot(
    df_ode.S,
    df_ode.I,
    label = "Фазовая траектория",
    xlabel = L"S(t)",
    ylabel = L"I(t)",
    title = "Фазовый портрет SIR модели",
    color = :blue,
    linewidth = 2,
    grid = true,
    size = (800, 500),
    legend = :topright,
)

# Добавим стрелки, показывающие направление движения системы
# по фазовой траектории.

for i in 1:50:length(df_ode.S)-1
    plot!(
        plt5,
        [df_ode.S[i], df_ode.S[i+1]],
        [df_ode.I[i], df_ode.I[i+1]],
        arrow = :closed,
        color = :blue,
        alpha = 0.5,
        label = false,
    )
end

# ## Эффективное репродуктивное число
#
# В процессе эпидемии число восприимчивых уменьшается.
# Поэтому используется эффективное репродуктивное число:
#
# $$
# R_e(t) = R_0\frac{S(t)}{N}.
# $$

df_ode[!, :Re] =
    R0 .* df_ode.S ./ df_ode.N

plt6 = @df df_ode plot(
    :t,
    :Re,
    label = L"R_e(t)",
    xlabel = "Время, дни",
    ylabel = L"R_e",
    title = "Динамика эффективного репродуктивного числа",
    color = :green,
    linewidth = 2,
    grid = true,
    size = (800, 400),
)

# Значение `Rₑ = 1` является пороговым.

hline!(
    plt6,
    [1.0],
    color = :red,
    linestyle = :dash,
    label = "Порог эпидемии (Rₑ=1)",
    linewidth = 1.5,
)

# Найдём момент, когда эффективное репродуктивное число
# впервые становится меньше единицы.

cross_idx = findfirst(x -> x < 1, df_ode.Re)

if !isnothing(cross_idx) && cross_idx > 1
    cross_time = df_ode.t[cross_idx]

    vline!(
        plt6,
        [cross_time],
        color = :black,
        linestyle = :dash,
        label = false,
        linewidth = 1,
    )

    annotate!(
        plt6,
        cross_time,
        1.2,
        text(
            "Rₑ<1 с $(round(cross_time, digits=1)) дня",
            8,
            :left,
        ),
    )
end

# ## Сводная панель результатов
#
# Для удобного сравнения основных характеристик модели
# разместим несколько графиков на одной панели.

plt7 = plot(
    layout = (2, 3),
    size = (1200, 800),
)

# ### Верхний ряд

plot!(
    plt7[1],
    df_ode.t,
    df_ode.S,
    label = L"S(t)",
    color = 1,
    linewidth = 2,
    title = "Восприимчивые",
)

plot!(
    plt7[2],
    df_ode.t,
    df_ode.I,
    label = L"I(t)",
    color = 2,
    linewidth = 2,
    title = "Зараженные",
)

plot!(
    plt7[3],
    df_ode.t,
    df_ode.R,
    label = L"R(t)",
    color = 3,
    linewidth = 2,
    title = "Выздоровевшие",
)

# ### Нижний ряд

plot!(
    plt7[4],
    df_ode.t,
    df_ode.I,
    label = L"I(t)",
    color = 2,
    linewidth = 2,
    yscale = :log10,
    title = "Лог. масштаб",
)

plot!(
    plt7[5],
    df_ode.S,
    df_ode.I,
    label = false,
    color = 4,
    linewidth = 2,
    title = "Фазовый портрет",
    xlabel = L"S",
    ylabel = L"I",
)

plot!(
    plt7[6],
    df_ode.t,
    df_ode.Re,
    label = L"R_e",
    color = :green,
    linewidth = 2,
    title = L"R_e(t)",
    hline = [1.0],
    linestyle = :dash,
    linecolor = :red,
)

# ## Сохранение графиков
#
# Сохраним все построенные графики в каталог
# `plots/sir_ode/`.

savefig(plt1, plotsdir(script_name, "sir_main.png"))
savefig(plt2, plotsdir(script_name, "sir_infected.png"))
savefig(plt3, plotsdir(script_name, "sir_log_scale.png"))
savefig(plt4, plotsdir(script_name, "sir_percentages.png"))
savefig(plt5, plotsdir(script_name, "sir_phase_portrait.png"))
savefig(plt6, plotsdir(script_name, "sir_effective_R.png"))
savefig(plt7, plotsdir(script_name, "sir_panel.png"))

# ## Оценка производительности
#
# Используем `BenchmarkTools` для измерения производительности
# численного решения системы.

println("\nБенчмарк решения:")

@benchmark solve(prob_ode, dt = δt)

# ## Анализ результатов
#
# Выведем основные количественные характеристики эпидемического процесса.

println("\n=== АНАЛИЗ РЕЗУЛЬТАТОВ ===")

println(
    "Общая численность популяции (контроль): N = ",
    round(df_ode.N[1], digits = 1),
)

println(
    "Пиковое число зараженных: I_max = ",
    round(peak_value, digits = 1),
)

println(
    "Время достижения пика: t_peak = ",
    round(peak_time, digits = 1),
    " дней",
)

println(
    "Итоговое число переболевших: R(∞) = ",
    round(df_ode.R[end], digits = 1),
)

println(
    "Доля переболевших: ",
    round(
        df_ode.R[end] / df_ode.N[1] * 100,
        digits = 1,
    ),
    "%",
)

# Для `R₀ > 1` дополнительно выведем теоретические характеристики модели.

if R0 > 1
    println("\nТеоретический анализ:")

    println(
        " - Порог коллективного иммунитета: ",
        round((1 - 1 / R0) * 100, digits = 1),
        "%",
    )

    println(
        " - Теоретический пик при S/N = 1/R0 = ",
        round(1 / R0, digits = 3),
    )
end
