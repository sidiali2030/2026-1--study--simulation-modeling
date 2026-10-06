##!/usr/bin/env julia
## add_packages.jl
using Pkg
Pkg.activate(".") # Активируем текущий проект
## ОСНОВНЫЕ ПАКЕТЫ ДЛЯ РАБОТЫ
packages = [
    # Агентное моделирование
    "Agents",
    
    # Общие пакеты
    "Distributions",
    
    "StatsPlots",
    
    "Dates",

    # Многокритериальная оптимизация параметров модели
    "BlackBoxOptim",

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
