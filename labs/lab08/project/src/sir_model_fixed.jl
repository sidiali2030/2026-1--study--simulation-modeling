using ResumableFunctions
using ConcurrentSim
using Distributions
using DataFrames
using Random


function increment_fixed!(a::Array{Int64})
    push!(a, a[length(a)] + 1)
end

function decrement_fixed!(a::Array{Int64})
    push!(a, a[length(a)] - 1)
end

function carryover_fixed!(a::Array{Int64})
    push!(a, a[length(a)])
end


mutable struct SIRPersonFixed
    id::Int64
    status::Symbol
end


mutable struct SIRModelFixed
    sim::ConcurrentSim.Simulation
    β::Float64
    c::Float64
    γ::Float64

    ta::Array{Float64}
    Sa::Array{Int64}
    Ia::Array{Int64}
    Ra::Array{Int64}

    allIndividuals::Array{SIRPersonFixed}
end


function infection_update_fixed!(
    sim::ConcurrentSim.Simulation,
    m::SIRModelFixed
)
    push!(m.ta, ConcurrentSim.now(sim))

    decrement_fixed!(m.Sa)
    increment_fixed!(m.Ia)
    carryover_fixed!(m.Ra)
end


function recovery_update_fixed!(
    sim::ConcurrentSim.Simulation,
    m::SIRModelFixed
)
    push!(m.ta, ConcurrentSim.now(sim))

    carryover_fixed!(m.Sa)
    decrement_fixed!(m.Ia)
    increment_fixed!(m.Ra)
end


@resumable function live_fixed(
    env::ConcurrentSim.Simulation,
    individual::SIRPersonFixed,
    m::SIRModelFixed
)

    while individual.status == :S

        @yield timeout(
            env,
            rand(Exponential(1 / m.c))
        )

        alter = individual

        while alter == individual

            N = length(m.allIndividuals)

            index = rand(
                DiscreteUniform(1, N)
            )

            alter = m.allIndividuals[index]
        end

        if alter.status == :I

            if rand(Uniform(0, 1)) < m.β

                individual.status = :I

                infection_update_fixed!(env, m)
            end
        end
    end


    if individual.status == :I

        # Детерминированная длительность болезни
        @yield timeout(
            env,
            1 / m.γ
        )

        individual.status = :R

        recovery_update_fixed!(env, m)
    end
end


function MakeSIRModelFixed(u0, p)

    (S, I, R) = u0

    N = S + I + R

    (β, c, γ) = p

    sim = ConcurrentSim.Simulation()

    allIndividuals = SIRPersonFixed[]

    for i = 1:S
        push!(
            allIndividuals,
            SIRPersonFixed(i, :S)
        )
    end

    for i = (S + 1):(S + I)
        push!(
            allIndividuals,
            SIRPersonFixed(i, :I)
        )
    end

    for i = (S + I + 1):N
        push!(
            allIndividuals,
            SIRPersonFixed(i, :R)
        )
    end

    ta = Float64[0.0]

    Sa = Int64[S]
    Ia = Int64[I]
    Ra = Int64[R]

    return SIRModelFixed(
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


function activate_fixed(m::SIRModelFixed)

    [
        @process live_fixed(
            m.sim,
            individual,
            m
        )
        for individual in m.allIndividuals
    ]
end


function sir_run_fixed(
    m::SIRModelFixed,
    tf::Float64
)
    ConcurrentSim.run(m.sim, tf)
end


function out_fixed(m::SIRModelFixed)

    result = DataFrame()

    result[!, :t] = m.ta
    result[!, :S] = m.Sa
    result[!, :I] = m.Ia
    result[!, :R] = m.Ra

    return result
end
