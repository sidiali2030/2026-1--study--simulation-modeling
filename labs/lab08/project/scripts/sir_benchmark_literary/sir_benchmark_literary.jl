using DrWatson
@quickactivate "project"

include(srcdir("sir_model.jl"))

using Random
using BenchmarkTools

tmax = 40.0

u0 = [
    9900,
    100,
    0
]

p = [
    0.05,
    10.0,
    0.25
]

function run_once()

    Random.seed!(1234)

    m = MakeSIRModel(
        u0,
        p
    )

    activate(m)

    sir_run(
        m,
        tmax
    )

    return nothing
end

println(
    "Benchmark для N = ",
    sum(u0)
)

result = @benchmark run_once()

display(result)
