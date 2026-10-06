using ResumableFunctions
using ConcurrentSim
using Distributions
using DataFrames
using Random


# ============================================================
# Вспомогательные функции
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
# Структура индивида
# ============================================================

mutable struct SEIRPerson
    id::Int64

    # Возможные состояния:
    # :S - восприимчивый
    # :E - инфицированный, но ещё не заразный
    # :I - инфекционный
    # :R - выздоровевший
    status::Symbol
end


# ============================================================
# Структура модели SEIR
# ============================================================

mutable struct SEIRModel

    sim::ConcurrentSim.Simulation

    # Вероятность заражения при контакте
    β::Float64

    # Частота контактов
    c::Float64

    # Интенсивность перехода E -> I
    σ::Float64

    # Интенсивность выздоровления
    γ::Float64

    # Временные ряды
    ta::Array{Float64}

    Sa::Array{Int64}
    Ea::Array{Int64}
    Ia::Array{Int64}
    Ra::Array{Int64}

    # Все индивиды
    allIndividuals::Array{SEIRPerson}
end


# ============================================================
# S -> E
# Событие заражения
# ============================================================

function exposure_update!(
    sim::ConcurrentSim.Simulation,
    m::SEIRModel
)

    push!(
        m.ta,
        ConcurrentSim.now(sim)
    )

    # S уменьшается
    decrement!(m.Sa)

    # E увеличивается
    increment!(m.Ea)

    # I без изменений
    carryover!(m.Ia)

    # R без изменений
    carryover!(m.Ra)
end


# ============================================================
# E -> I
# Конец латентного периода
# ============================================================

function infectious_update!(
    sim::ConcurrentSim.Simulation,
    m::SEIRModel
)

    push!(
        m.ta,
        ConcurrentSim.now(sim)
    )

    # S без изменений
    carryover!(m.Sa)

    # E уменьшается
    decrement!(m.Ea)

    # I увеличивается
    increment!(m.Ia)

    # R без изменений
    carryover!(m.Ra)
end


# ============================================================
# I -> R
# Выздоровление
# ============================================================

function recovery_update!(
    sim::ConcurrentSim.Simulation,
    m::SEIRModel
)

    push!(
        m.ta,
        ConcurrentSim.now(sim)
    )

    # S без изменений
    carryover!(m.Sa)

    # E без изменений
    carryover!(m.Ea)

    # I уменьшается
    decrement!(m.Ia)

    # R увеличивается
    increment!(m.Ra)
end


# ============================================================
# Жизненный цикл индивида
# ============================================================

@resumable function live(
    env::ConcurrentSim.Simulation,
    individual::SEIRPerson,
    m::SEIRModel
)

    # --------------------------------------------------------
    # Состояние S
    # --------------------------------------------------------

    while individual.status == :S

        # Время до следующего контакта
        @yield timeout(
            env,
            rand(
                Exponential(1 / m.c)
            )
        )

        # Выбираем другого индивида
        alter = individual

        while alter == individual

            N = length(
                m.allIndividuals
            )

            index = rand(
                DiscreteUniform(1, N)
            )

            alter =
                m.allIndividuals[index]
        end


        # Заразить может только I
        if alter.status == :I

            if rand(
                Uniform(0, 1)
            ) < m.β

                # В SEIR после заражения:
                # S -> E
                individual.status = :E

                exposure_update!(
                    env,
                    m
                )
            end
        end
    end


    # --------------------------------------------------------
    # Состояние E
    # Латентный период
    # --------------------------------------------------------

    if individual.status == :E

        # По заданию 8.5.7:
        # ожидание экспоненциально распределено
        # со средним 1 / σ

        @yield timeout(
            env,
            rand(
                Exponential(1 / m.σ)
            )
        )

        # E -> I
        individual.status = :I

        infectious_update!(
            env,
            m
        )
    end


    # --------------------------------------------------------
    # Состояние I
    # --------------------------------------------------------

    if individual.status == :I

        # Время до выздоровления
        @yield timeout(
            env,
            rand(
                Exponential(1 / m.γ)
            )
        )

        # I -> R
        individual.status = :R

        recovery_update!(
            env,
            m
        )
    end
end


# ============================================================
# Создание SEIR-модели
#
# u0 = [S0, E0, I0, R0]
# p  = [β, c, σ, γ]
# ============================================================

function MakeSEIRModel(u0, p)

    (S, E, I, R) = u0

    N = S + E + I + R

    (β, c, σ, γ) = p

    sim =
        ConcurrentSim.Simulation()

    allIndividuals =
        SEIRPerson[]


    # --------------------------------------------------------
    # S
    # --------------------------------------------------------

    for i = 1:S

        push!(
            allIndividuals,
            SEIRPerson(
                i,
                :S
            )
        )
    end


    # --------------------------------------------------------
    # E
    # --------------------------------------------------------

    for i = (S + 1):(S + E)

        push!(
            allIndividuals,
            SEIRPerson(
                i,
                :E
            )
        )
    end


    # --------------------------------------------------------
    # I
    # --------------------------------------------------------

    for i =
        (S + E + 1):(S + E + I)

        push!(
            allIndividuals,
            SEIRPerson(
                i,
                :I
            )
        )
    end


    # --------------------------------------------------------
    # R
    # --------------------------------------------------------

    for i =
        (S + E + I + 1):N

        push!(
            allIndividuals,
            SEIRPerson(
                i,
                :R
            )
        )
    end


    # --------------------------------------------------------
    # Начальные значения статистики
    # --------------------------------------------------------

    ta = Float64[0.0]

    Sa = Int64[S]
    Ea = Int64[E]
    Ia = Int64[I]
    Ra = Int64[R]


    return SEIRModel(
        sim,
        β,
        c,
        σ,
        γ,
        ta,
        Sa,
        Ea,
        Ia,
        Ra,
        allIndividuals
    )
end


# ============================================================
# Активация всех процессов
# ============================================================

function activate(
    m::SEIRModel
)

    [
        @process live(
            m.sim,
            individual,
            m
        )

        for individual
        in m.allIndividuals
    ]
end


# ============================================================
# Запуск модели
# ============================================================

function seir_run(
    m::SEIRModel,
    tf::Float64
)

    ConcurrentSim.run(
        m.sim,
        tf
    )
end


# ============================================================
# Получение результатов
# ============================================================

function out(
    m::SEIRModel
)

    result =
        DataFrame()

    result[!, :t] =
        m.ta

    result[!, :S] =
        m.Sa

    result[!, :E] =
        m.Ea

    result[!, :I] =
        m.Ia

    result[!, :R] =
        m.Ra

    return result
end
