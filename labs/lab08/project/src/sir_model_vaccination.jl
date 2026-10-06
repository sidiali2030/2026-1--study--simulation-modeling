using ResumableFunctions
using ConcurrentSim
using Distributions
using DataFrames
using Random

# ============================================================
# Вспомогательные функции для обновления массивов состояния
# ============================================================

function increment!(a::Array{Int64})
    push!(a, a[end] + 1)
end

function decrement!(a::Array{Int64})
    push!(a, a[end] - 1)
end

function carryover!(a::Array{Int64})
    push!(a, a[end])
end


# ============================================================
# Структуры данных
# ============================================================

mutable struct SIRPerson
    id::Int64
    status::Symbol       # :S, :I, :R
end


mutable struct SIRModel
    sim::ConcurrentSim.Simulation

    β::Float64
    c::Float64
    γ::Float64

    ta::Array{Float64}
    Sa::Array{Int64}
    Ia::Array{Int64}
    Ra::Array{Int64}

    allIndividuals::Array{SIRPerson}
end


# ============================================================
# Обновление статистики при заражении
# ============================================================

function infection_update!(
    sim::ConcurrentSim.Simulation,
    m::SIRModel
)
    push!(m.ta, ConcurrentSim.now(sim))

    decrement!(m.Sa)
    increment!(m.Ia)
    carryover!(m.Ra)
end


# ============================================================
# Обновление статистики при выздоровлении
# ============================================================

function recovery_update!(
    sim::ConcurrentSim.Simulation,
    m::SIRModel
)
    push!(m.ta, ConcurrentSim.now(sim))

    carryover!(m.Sa)
    decrement!(m.Ia)
    increment!(m.Ra)
end


# ============================================================
# Обновление статистики при вакцинации
# S -> R
# ============================================================

function vaccination_update!(
    sim::ConcurrentSim.Simulation,
    m::SIRModel
)
    push!(m.ta, ConcurrentSim.now(sim))

    decrement!(m.Sa)
    carryover!(m.Ia)
    increment!(m.Ra)
end


# ============================================================
# Основной жизненный цикл индивида
# ============================================================

@resumable function live(
    env::ConcurrentSim.Simulation,
    individual::SIRPerson,
    m::SIRModel
)

    # Пока индивид восприимчив
    while individual.status == :S

        # Время до следующего контакта
        @yield timeout(
            env,
            rand(Exponential(1 / m.c))
        )

        # Важно:
        # во время ожидания человек мог быть вакцинирован.
        # В таком случае процесс заражения прекращаем.
        if individual.status != :S
            break
        end

        alter = individual

        # Выбираем другого человека
        while alter == individual

            N = length(m.allIndividuals)

            index = rand(
                DiscreteUniform(1, N)
            )

            alter = m.allIndividuals[index]
        end

        # Если контакт с инфицированным
        if alter.status == :I

            if rand(Uniform(0, 1)) < m.β

                individual.status = :I

                infection_update!(env, m)
            end
        end
    end


    # Если человек инфицирован
    if individual.status == :I

        # Стохастическая длительность болезни
        @yield timeout(
            env,
            rand(Exponential(1 / m.γ))
        )

        individual.status = :R

        recovery_update!(env, m)
    end
end


# ============================================================
# Процесс вакцинации
# ============================================================

@resumable function vaccinate(
    env::ConcurrentSim.Simulation,
    m::SIRModel,
    time::Float64,
    fraction::Float64
)

    # Ждём момента начала вакцинации
    @yield timeout(env, time)

    # Получаем всех восприимчивых людей
    susceptible = [
        individual
        for individual in m.allIndividuals
        if individual.status == :S
    ]

    S_current = length(susceptible)

    # Количество вакцинируемых
    n_vaccinated = round(
        Int,
        fraction * S_current
    )

    # Не позволяем выйти за пределы
    n_vaccinated = clamp(
        n_vaccinated,
        0,
        S_current
    )

    # Перемешиваем список, чтобы выбрать случайных людей
    shuffle!(susceptible)

    println()
    println("========================================")
    println("ВАКЦИНАЦИЯ")
    println("========================================")
    println("Время вакцинации: ", ConcurrentSim.now(env))
    println("Восприимчивых до вакцинации: ", S_current)
    println("Доля вакцинации: ", fraction)
    println("Количество вакцинируемых: ", n_vaccinated)

    # Переводим выбранных людей S -> R
    for i in 1:n_vaccinated

        individual = susceptible[i]

        # Дополнительная проверка
        if individual.status == :S

            individual.status = :R

            vaccination_update!(env, m)
        end
    end

    println("S после вакцинации: ", m.Sa[end])
    println("I после вакцинации: ", m.Ia[end])
    println("R после вакцинации: ", m.Ra[end])
    println("========================================")
end


# ============================================================
# Создание модели
# ============================================================

function MakeSIRModel(u0, p)

    (S, I, R) = u0

    N = S + I + R

    (β, c, γ) = p

    sim = ConcurrentSim.Simulation()

    allIndividuals = SIRPerson[]


    # Восприимчивые
    for i = 1:S

        push!(
            allIndividuals,
            SIRPerson(i, :S)
        )
    end


    # Инфицированные
    for i = (S + 1):(S + I)

        push!(
            allIndividuals,
            SIRPerson(i, :I)
        )
    end


    # Выздоровевшие
    for i = (S + I + 1):N

        push!(
            allIndividuals,
            SIRPerson(i, :R)
        )
    end


    ta = Float64[0.0]

    Sa = Int64[S]
    Ia = Int64[I]
    Ra = Int64[R]


    return SIRModel(
        sim,
        β,
        c,
        γ,
        ta,
        Sa,
        Ia,
        Ra,
        allIndividuals
    )
end


# ============================================================
# Активация процессов индивидов
# ============================================================

function activate(m::SIRModel)

    [
        @process live(
            m.sim,
            individual,
            m
        )
        for individual in m.allIndividuals
    ]
end


# ============================================================
# Активация вакцинации
# ============================================================

function activate_vaccination!(
    m::SIRModel,
    time::Float64,
    fraction::Float64
)

    @process vaccinate(
        m.sim,
        m,
        time,
        fraction
    )
end


# ============================================================
# Запуск симуляции
# ============================================================

function sir_run(
    m::SIRModel,
    tf::Float64
)

    ConcurrentSim.run(
        m.sim,
        tf
    )
end


# ============================================================
# Формирование результата
# ============================================================

function out(m::SIRModel)

    result = DataFrame()

    result[!, :t] = m.ta
    result[!, :S] = m.Sa
    result[!, :I] = m.Ia
    result[!, :R] = m.Ra

    return result
end
