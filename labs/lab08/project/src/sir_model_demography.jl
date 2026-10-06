using ResumableFunctions
using ConcurrentSim
using Distributions
using DataFrames
using Random


# Вспомогательные функции
function increment!(a::Array{Int64})
    push!(a, a[length(a)] + 1)
end

function decrement!(a::Array{Int64})
    push!(a, a[length(a)] - 1)
end

function carryover!(a::Array{Int64})
    push!(a, a[length(a)])
end


# Индивид
mutable struct SIRPerson
    id::Int64
    status::Symbol # :S, :I, :R, :D
end


# Модель SIR с демографическими событиями
mutable struct SIRModel
    sim::ConcurrentSim.Simulation

    β::Float64
    c::Float64
    γ::Float64

    μ::Float64
    λ_birth::Float64

    ta::Array{Float64}
    Sa::Array{Int64}
    Ia::Array{Int64}
    Ra::Array{Int64}
    Da::Array{Int64}

    allIndividuals::Array{SIRPerson}
end


# Заражение
function infection_update!(
    sim::ConcurrentSim.Simulation,
    m::SIRModel
)
    push!(m.ta, ConcurrentSim.now(sim))

    decrement!(m.Sa)
    increment!(m.Ia)
    carryover!(m.Ra)
    carryover!(m.Da)
end


# Выздоровление
function recovery_update!(
    sim::ConcurrentSim.Simulation,
    m::SIRModel
)
    push!(m.ta, ConcurrentSim.now(sim))

    carryover!(m.Sa)
    decrement!(m.Ia)
    increment!(m.Ra)
    carryover!(m.Da)
end


# Смерть
function death_update!(
    sim::ConcurrentSim.Simulation,
    m::SIRModel,
    old_status::Symbol
)
    push!(m.ta, ConcurrentSim.now(sim))

    if old_status == :S

        decrement!(m.Sa)
        carryover!(m.Ia)
        carryover!(m.Ra)

    elseif old_status == :I

        carryover!(m.Sa)
        decrement!(m.Ia)
        carryover!(m.Ra)

    elseif old_status == :R

        carryover!(m.Sa)
        carryover!(m.Ia)
        decrement!(m.Ra)
    end

    increment!(m.Da)
end


# Рождение
function birth_update!(
    sim::ConcurrentSim.Simulation,
    m::SIRModel
)
    push!(m.ta, ConcurrentSim.now(sim))

    increment!(m.Sa)
    carryover!(m.Ia)
    carryover!(m.Ra)
    carryover!(m.Da)
end


# Жизненный цикл: заражение и выздоровление
@resumable function live(
    env::ConcurrentSim.Simulation,
    individual::SIRPerson,
    m::SIRModel
)
    while individual.status == :S

        @yield timeout(
            env,
            rand(Exponential(1 / m.c))
        )

        if individual.status != :S
            break
        end

        livingIndividuals = [
            person for person in m.allIndividuals
            if person.status != :D
        ]

        if length(livingIndividuals) <= 1
            continue
        end

        alter = individual

        while alter == individual

            N = length(livingIndividuals)

            index = rand(
                DiscreteUniform(1, N)
            )

            alter = livingIndividuals[index]
        end

        if alter.status == :I

            if rand(Uniform(0, 1)) < m.β

                individual.status = :I

                infection_update!(
                    env,
                    m
                )
            end
        end
    end

    if individual.status == :I

        @yield timeout(
            env,
            rand(Exponential(1 / m.γ))
        )

        if individual.status == :I

            individual.status = :R

            recovery_update!(
                env,
                m
            )
        end
    end
end


# Процесс смерти
@resumable function die(
    env::ConcurrentSim.Simulation,
    individual::SIRPerson,
    m::SIRModel
)
    @yield timeout(
        env,
        rand(Exponential(1 / m.μ))
    )

    if individual.status != :D

        old_status = individual.status

        individual.status = :D

        death_update!(
            env,
            m,
            old_status
        )
    end
end


# Процесс рождения
@resumable function births(
    env::ConcurrentSim.Simulation,
    m::SIRModel
)
    while true

        @yield timeout(
            env,
            rand(Exponential(1 / m.λ_birth))
        )

        new_id = length(m.allIndividuals) + 1

        new_person = SIRPerson(
            new_id,
            :S
        )

        push!(
            m.allIndividuals,
            new_person
        )

        birth_update!(
            env,
            m
        )

        @process live(
            env,
            new_person,
            m
        )

        @process die(
            env,
            new_person,
            m
        )
    end
end


# Создание модели
function MakeSIRModel(u0, p)

    (S, I, R) = u0
    N = S + I + R

    (β, c, γ, μ, λ_birth) = p

    sim = ConcurrentSim.Simulation()

    allIndividuals = SIRPerson[]

    for i = 1:S
        push!(
            allIndividuals,
            SIRPerson(i, :S)
        )
    end

    for i = (S + 1):(S + I)
        push!(
            allIndividuals,
            SIRPerson(i, :I)
        )
    end

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
    Da = Int64[0]

    SIRModel(
        sim,
        β,
        c,
        γ,
        μ,
        λ_birth,
        ta,
        Sa,
        Ia,
        Ra,
        Da,
        allIndividuals
    )
end


# Активация процессов
function activate(m::SIRModel)

    for individual in m.allIndividuals

        @process live(
            m.sim,
            individual,
            m
        )

        @process die(
            m.sim,
            individual,
            m
        )
    end

    @process births(
        m.sim,
        m
    )
end


# Запуск
function sir_run(
    m::SIRModel,
    tf::Float64
)
    ConcurrentSim.run(
        m.sim,
        tf
    )
end


# Результаты
function out(m::SIRModel)

    result = DataFrame()

    result[!, :t] = m.ta
    result[!, :S] = m.Sa
    result[!, :I] = m.Ia
    result[!, :R] = m.Ra
    result[!, :D] = m.Da

    return result
end
