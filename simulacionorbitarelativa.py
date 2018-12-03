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
for x in [row[1] for row in trayectorias]:
    for y in [row[3] for row in trayectorias]:
        if y - x > max_x:
            max_x = y - x
        elif y - x < min_x:
            min_x = y - x
for x in [row[2] for row in trayectorias]:
    for y in [row[4] for row in trayectorias]:
        if y - x > max_y:
            max_y = y - x
        elif y - x < min_y:
            min_y = y - x

print(min_x, max_x, min_y, max_y)

fig = plt.figure()
ax = plt.axes(xlim=(min_x, max_x), ylim=(min_y, max_y))
print(radios[0] / 10 ** (np.log10(radios[0]) - 1))
sun, = ax.plot([0, 0], [0, 0], '-o', markersize=2 * radios[0] / 10 ** (np.log10(radios[0]) - 1))
earth, = ax.plot([], [], '-o', markersize=2 * radios[1] / 10 ** (np.log10(radios[0]) - 1))
earthtr, = ax.plot([], [], '-')


def init():
    earth.set_data([[trayectorias[0][3] - trayectorias[0][1]], [[trayectorias[0][4] - trayectorias[0][2]]]])
    earthtr.set_data([], [])
    return sun, earth


def animate(i):
    earthtr.set_data([(row[3] - row[1]) for row in trayectorias[0:i + 1]], [(row[4] - row[2]) for row in trayectorias[0:i + 1]])
    sunx = trayectorias[i + 1][1]
    suny = trayectorias[i + 1][2]
    earthx = trayectorias[i + 1][3]
    earthy = trayectorias[i + 1][4]
    earth.set_data([earthx - sunx], [earthy - suny])
    return earth, earthtr


anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=250 * 3, interval=15, blit=True)

plt.show()
