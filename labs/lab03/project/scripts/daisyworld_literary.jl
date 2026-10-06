# # Базовая визуализация модели Daisyworld
#
# В данном скрипте выполняется базовая визуализация агентной модели
# Daisyworld. Модель строится на основе ранее реализованного файла
# `src/daisyworld.jl`.
#
# Цель визуализации — показать распределение температуры поверхности
# в виде тепловой карты и отобразить чёрные и белые маргаритки
# в соответствии с их типом.

using DrWatson
@quickactivate "project"
using Agents
using DataFrames
using Plots
include(srcdir("daisyworld.jl"))
using CairoMakie

# ## Инициализация модели
#
# Создаём экземпляр модели Daisyworld с параметрами по умолчанию.

model = daisyworld()

# ## Настройка отображения агентов
#
# Цвет каждого агента определяется его типом `breed`.
# Чёрные и белые маргаритки будут отображаться в соответствии
# со значением этого свойства.

daisycolor(a::Daisy) = a.breed

# ## Параметры визуализации
#
# Задаём общие параметры для построения графиков.
#
# - `agent_color` определяет цвет маргаритки;
# - `agent_size` задаёт размер символа агента;
# - `agent_marker` задаёт символ маргаритки;
# - `heatarray` указывает, что тепловая карта строится по температуре;
# - `heatkwargs` задаёт диапазон температур от -20 до 60.

plotkwargs = (
agent_color=daisycolor, agent_size = 20, agent_marker = '✿',
heatarray = :temperature,
heatkwargs = (colorrange = (-20, 60),),
)

# ## Начальное состояние модели
#
# Строим визуализацию модели до выполнения шагов симуляции.

plt1, _ = abmplot(model; plotkwargs...)

# ## Состояние после 5 шагов
#
# Продвигаем модель на 5 шагов и строим следующую визуализацию.

step!(model, 5)
plt2, _ = abmplot(model; heatarray = model.temperature, plotkwargs...)

# ## Состояние после следующих 40 шагов
#
# Выполняем ещё 40 шагов модели и строим третью визуализацию.

step!(model, 40)
plt3, _ = abmplot(model; heatarray = model.temperature, plotkwargs...)

# ## Сохранение результатов
#
# Полученные изображения сохраняются в каталог `plots`,
# который определяется структурой проекта DrWatson.

save(plotsdir("daisy_step001.png"), plt1)
save(plotsdir("daisy_step005.png"), plt2)
save(plotsdir("daisy_step040.png"), plt3)
