##!/usr/bin/env julia
## add_packages.jl
using Pkg
Pkg.activate(".") # Активируем текущий проект
## ОСНОВНЫЕ ПАКЕТЫ ДЛЯ РАБОТЫ
packages = [
    # Дискретно-событийное моделирование
    "ResumableFunctions",
    "ConcurrentSim",


    # Агентное моделирование
    "Agents",
    "CairoMakie",

    # Общие пакеты
    "Distributions",
    "StatsPlots",
    "Random",
    "Dates",

    # Организация проекта
    "DrWatson",

    # Дифференциальные уравнения
    "DifferentialEquations",

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
