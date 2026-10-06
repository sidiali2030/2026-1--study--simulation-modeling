##!/usr/bin/env julia
## add_packages.jl
using Pkg
Pkg.activate(".") # Активируем текущий проект
## ОСНОВНЫЕ ПАКЕТЫ ДЛЯ РАБОТЫ
packages = [
    # Общие пакеты
    "Distributions",
    "StatsPlots",
    "StableRNGs",
    "StatsBase",
    "Dates",
    "ConcurrentSim",
    "ResumableFunctions",

    # Организация проекта
    "DrWatson",

    # Визуализация
    "Plots",

    # Данные
    "DataFrames",
    "CSV",
    "JLD2",

    # Literate programming
    "Literate",
    "IJulia",

    # Тестирование производительности
    "BenchmarkTools",
    
    # Отчёты
    "Quarto"
]

println("Установка базовых пакетов...")
Pkg.add(packages)
println("\n Все пакеты установлены!")
println("Для проверки: using DrWatson, DifferentialEquations, Plots")
