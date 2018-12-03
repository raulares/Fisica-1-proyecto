# coding: utf-8

from matplotlib import pyplot as plt
from matplotlib import animation
import numpy as np

datos = open('casosimple1.dump', 'r')
tabla = datos.read().splitlines()
datos.close

n = int(tabla[1])
paso = int(tabla[4])
planetas = tabla[7:7 + n]
masas = [float(x) for x in tabla[9 + n:9 + 2 * n]]
radios = [float(x) for x in tabla[11 + 2 * n:11 + 3 * n]]
trayectorias = [x.split() for x in tabla[13 + 3 * n:len(tabla)]]
trayectorias = [[float(y) for y in x] for x in trayectorias]

max_x = 0
min_x = 0
max_y = 0
min_y = 0
for x in [row[1] for row in trayectorias] + [row[3] for row in trayectorias]:
    if x > max_x:
        max_x = x
    elif x < min_x:
        min_x = x
for y in [row[2] for row in trayectorias] + [row[4] for row in trayectorias]:
    if y > max_y:
        max_y = y
    elif y < min_y:
        min_y = y

fig = plt.figure()
ax = plt.axes(xlim=(min_x, max_x), ylim=(min_y, max_y))
print(radios[0] / 10 ** (np.log10(radios[0]) - 1))
sun, = ax.plot([], [], '-o', markersize=2 * radios[0] / 10 ** (np.log10(radios[0]) - 1))
suntr, = ax.plot([], [], '-')
earth, = ax.plot([], [], '-o', markersize=2 * radios[1] / 10 ** (np.log10(radios[0]) - 1))
earthtr, = ax.plot([], [], '-')


def init():
    sun.set_data([trayectorias[0][1]], [[trayectorias[0][2]]])
    earth.set_data([[trayectorias[0][3]], [[trayectorias[0][4]]]])
    suntr.set_data([], [])
    earthtr.set_data([], [])
    return sun, earth


def animate(i):
    suntr.set_data([row[1] for row in trayectorias[0:i + 1]], [row[2] for row in trayectorias[0:i + 1]])
    earthtr.set_data([row[3] for row in trayectorias[0:i + 1]], [row[4] for row in trayectorias[0:i + 1]])
    sunx = [trayectorias[i + 1][1]]
    suny = [trayectorias[i + 1][2]]
    earthx = [trayectorias[i + 1][3]]
    earthy = [trayectorias[i + 1][4]]
    sun.set_data(sunx, suny)
    earth.set_data(earthx, earthy)
    return sun, earth, suntr, earthtr


anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=250 * 5, interval=10, blit=True)

plt.show()
