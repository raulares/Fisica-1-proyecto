# coding: utf-8

import matplotlib


class PlotTraj:

    def __init__(self, filename):
        archivo = open(filename, 'r')
        datos = archivo.read().splitlines()

        self.cuerpos = int(datos[1])
        self.paso = int(datos[4])
        self.planetas = datos[7:7 + self.cuerpos]
        self.masas = [float(x) for x in datos[9 + self.cuerpos:9 + 2 * self.cuerpos]]
        self.radios = [float(x) for x in datos[11 + 2 * self.cuerpos:11 + 3 * self.cuerpos]]
        self.trayectorias = [x.split() for x in datos[13 + 3 * self.cuerpos:len(datos)]]
        self.trayectorias = [[float(y) for y in x] for x in self.trayectorias]
        archivo.close

    def _update(self, frame):


SolTierra = PlotTraj('sun_earth.dump')

print(SolTierra.trayectorias)
print(len(SolTierra.trayectorias))
