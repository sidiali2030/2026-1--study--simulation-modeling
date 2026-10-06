##!/usr/bin/env julia
## add_packages.jl
using Pkg
Pkg.activate(".") # Активируем текущий проект
## ОСНОВНЫЕ ПАКЕТЫ ДЛЯ РАБОТЫ
packages = [
    # Агентное моделирование
    "Agents",
    "AlgebraicPetri",
    
    # Категориальная алгебра и визуализация сетей Петри
    "Catlab",
    
    # Общие пакеты
    "Distributions",
    "StatsPlots",
    "StatsBase",
    "Dates",

    # Многокритериальная оптимизация параметров модели
    "BlackBoxOptim",

    # Организация проекта
    "DrWatson",

    # Дифференциальные уравнения
    "DifferentialEquations",
    "OrdinaryDiffEq",

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
