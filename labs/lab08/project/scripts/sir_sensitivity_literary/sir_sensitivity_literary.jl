using DrWatson
@quickactivate "project"

include(srcdir("sir_model.jl"))

using Random
using DataFrames
using StatsPlots

tmax = 40.0

u0 = [
    990,
    10,
    0
]

betas = [
    0.03,
    0.05,
    0.07
]

for β in betas


    Random.seed!(1234)

    p = [
        β,
        10.0,
        0.25
    ]


    m = MakeSIRModel(
        u0,
        p
    )


    activate(m)


    sir_run(
        m,
        tmax
    )


    data = out(m)


    peak_I = maximum(
        data.I
    )

    peak_index = argmax(
        data.I
    )


    peak_time =
        data.t[peak_index]


    final_R =
        data.R[end]

    println(
        "β = $β, ",
        "peak I = $peak_I, ",
        "peak time = $peak_time, ",
        "final R = $final_R"
    )

    pfig = @df data plot(
        :t,
        [:S :I :R],
        label = ["S" "I" "R"],
        xlabel = "Время",
        ylabel = "Численность",
        title = "β = $β"
    )

    savefig(
        pfig,
        plotsdir(
            "sir_beta_$β.png"
        )
    )
end
