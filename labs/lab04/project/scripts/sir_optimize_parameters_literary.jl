# # Многокритериальная оптимизация параметров модели SIR
#
# В данном эксперименте выполняется многокритериальная оптимизация
# параметров агентной модели SIR.
#
# Цель оптимизации состоит в одновременном уменьшении двух показателей:
#
# - максимальной доли инфицированных агентов;
# - доли умерших агентов.
#
# В качестве изменяемых параметров рассматриваются:
#
# - коэффициент передачи инфекции `β_und`;
# - время выявления инфекции `detection_time`;
# - вероятность смерти `death_rate`.
#
# Для поиска компромиссных решений используется многокритериальный
# эволюционный алгоритм Borg MOEA из пакета `BlackBoxOptim`.

# ## Подключение библиотек
#
# Активируем окружение проекта DrWatson и подключаем библиотеки,
# необходимые для оптимизации, статистической обработки результатов
# и работы с агентной моделью.

using DrWatson
@quickactivate "project"

using BlackBoxOptim
using Random
using Statistics
using Agents
using JLD2

# Подключаем реализацию агентной модели SIR.

include(srcdir("sir_model.jl"))

# ## Целевая функция многокритериальной оптимизации
#
# Определим функцию `cost_multi`, которая принимает вектор
# оптимизируемых параметров `x`.
#
# Компоненты вектора интерпретируются следующим образом:
#
# - `x[1]` — коэффициент передачи инфекции `β_und`;
# - `x[2]` — время выявления инфекции `detection_time`;
# - `x[3]` — вероятность смерти `death_rate`.
#
# Функция возвращает два значения:
#
# 1. среднюю пиковую долю инфицированных;
# 2. среднюю долю умерших.
#
# Оба показателя должны быть минимизированы.

function cost_multi(x)

    

    β_und = fill(x[1], 3)
    β_det = fill(x[1] / 10, 3)

   

    detection_time = round(Int, x[2])

   

    death_rate = x[3]

    
    infected_frac(model) =
        count(a.status == :I for a in allagents(model)) / nagents(model)

   

    dead_count(model) =
        3000 - nagents(model)

    

    replicates = 5

    peak_vals = Float64[]
    dead_vals = Int[]

    

    for rep = 1:replicates

        

        model = initialize_sir(
            ;
            Ns = [1000, 1000, 1000],
            β_und = β_und,
            β_det = β_det,
            infection_period = 14,
            detection_time = detection_time,
            death_rate = death_rate,
            reinfection_probability = 0.1,
            Is = [0, 0, 1],
            seed = 42 + rep,
            n_steps = 100,
        )

        

        peak_infected = 0.0

        

        for step = 1:100

            Agents.step!(model, 1)

            frac = infected_frac(model)

            if frac > peak_infected
                peak_infected = frac
            end
        end

        

        push!(peak_vals, peak_infected)
        push!(dead_vals, dead_count(model))
    end

   
    mean_peak = mean(peak_vals)
    mean_death_fraction = mean(dead_vals) / 3000

    

    return (
        mean_peak,
        mean_death_fraction,
    )
end

# ## Диапазоны оптимизируемых параметров
#
# Оптимизация выполняется в следующих пределах:
#
# - `β_und` от `0.1` до `1.0`;
# - `detection_time` от `3` до `14` дней;
# - `death_rate` от `0.01` до `0.1`.

search_range = [
    (0.1, 1.0),
    (3.0, 14.0),
    (0.01, 0.1),
]

# ## Запуск многокритериальной оптимизации
#
# Используем алгоритм `borg_moea`.
#
# Поскольку одновременно оптимизируются два показателя,
# задаётся схема `ParetoFitnessScheme{2}`.
#
# Значение `is_minimizing = true` указывает, что оба критерия
# необходимо минимизировать.
#
# Максимальное время работы алгоритма составляет 120 секунд.

result = bboptimize(
    cost_multi,
    Method = :borg_moea,
    FitnessScheme = ParetoFitnessScheme{2}(
        is_minimizing = true,
    ),
    SearchRange = search_range,
    NumDimensions = 3,
    MaxTime = 120,
    TraceMode = :compact,
)

# ## Получение результата оптимизации
#
# Из результата оптимизации извлекаем найденный набор параметров
# и соответствующие значения целевых функций.

best = best_candidate(result)
fitness = best_fitness(result)

# ## Вывод оптимальных параметров
#
# Выведем найденные значения параметров модели.

println("Оптимальные параметры:")
println("β_und = $(best[1])")
println(
    "Время выявления = $(round(Int, best[2])) дней",
)
println("Смертность = $(best[3])")

# Выведем достигнутые значения критериев.

println("Достигнутые показатели:")
println("Пик заболеваемости: $(fitness[1])")
println("Доля умерших: $(fitness[2])")

# ## Сохранение результатов
#
# Результат оптимизации сохраняется в файл JLD2.
#
# В файл записываются:
#
# - `best` — найденный набор параметров;
# - `fitness` — соответствующие значения целевых функций.

save(
    datadir("optimization_result.jld2"),
    Dict(
        "best" => best,
        "fitness" => fitness,
    ),
)

# ## Вывод
#
# В результате выполнена многокритериальная оптимизация параметров
# агентной модели SIR. Одновременно минимизировались пиковая
# заболеваемость и доля умерших.
#
# Использование нескольких независимых прогонов для каждого
# набора параметров позволяет уменьшить влияние случайности,
# характерной для агентного моделирования.
#
# Алгоритм Borg MOEA позволяет искать компромиссные решения
# между двумя конкурирующими критериями и определить набор
# параметров, обеспечивающий более благоприятную динамику
# эпидемического процесса.
