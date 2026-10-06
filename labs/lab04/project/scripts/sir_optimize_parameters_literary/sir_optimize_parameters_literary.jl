using DrWatson
@quickactivate "project"

using BlackBoxOptim
using Random
using Statistics
using Agents
using JLD2

include(srcdir("sir_model.jl"))

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

search_range = [
    (0.1, 1.0),
    (3.0, 14.0),
    (0.01, 0.1),
]

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

best = best_candidate(result)
fitness = best_fitness(result)

println("Оптимальные параметры:")
println("β_und = $(best[1])")
println(
    "Время выявления = $(round(Int, best[2])) дней",
)
println("Смертность = $(best[3])")

println("Достигнутые показатели:")
println("Пик заболеваемости: $(fitness[1])")
println("Доля умерших: $(fitness[2])")

save(
    datadir("optimization_result.jld2"),
    Dict(
        "best" => best,
        "fitness" => fitness,
    ),
)
