#include <iostream>
#include <cmath>
#include <string>

#define G 6.67408e-11
#define SOL_MASA 1.989e30
#define SOL_RADIO 695700000.0
#define TIERRA_X0 -147095000000.0
#define TIERRA_Y0 0.0
#define TIERRA_VX0 0.0
#define TIERRA_VY0 -30300.0
#define TIERRA_MASA 5.972e24
#define TIERRA_RADIO 6371000.0
#define PLUTON_MASA 1.25e22
#define PLUTON_RADIO 1185000.0
#define CARONTE_X0 -17535294.0
#define CARONTE_Y0 0.0
#define CARONTE_VX0 0.0
#define CARONTE_VY0 -230.0
#define CARONTE_MASA 1.52e21
#define CARONTE_RADIO 604000.0
#define DT 1800.0
#define NUM_CUERPOS 2
#define NUM_PASOS 365

using namespace std;

class Cuerpo {
    public:
        Cuerpo (string nombre, double x, double y, double vx, double vy, double m, double r);
        ~Cuerpo () {}

        void asignarNombre (string nombre) { suNombre = nombre; }
        void asignarPosX (double x) { suX = x; }
        void asignarPosY (double y) { suY = y; }
        void asignarVelX (double x) { suVX = x; }
        void asignarVelY (double y) { suVY = y; }
        void asignarMasa (double m) { suMasa = m; }
        void asignarRadio (double r) { suRadio = r; }

        string obtenerNombre () { return suNombre; }
        double x () { return suX; }
        double y () { return suY; }
        double vx () { return suVX; }
        double vy () { return suVY; }
        double obtenerMasa () { return suMasa; }
        double obtenerRadio () { return suRadio; }

        void calcular_aceleracion(double x, double y, double & ax, double & ay, Cuerpo * cuerpo2);
        void paso(double dt, Cuerpo * cuerpo);

    private:
        string suNombre;
        double suX;
        double suY;
        double suVX;
        double suVY;
        double suMasa;
        double suRadio;
};

Cuerpo::Cuerpo (string nombre, double x, double y, double vx, double vy, double m, double r) {
    suNombre = nombre;
    suX = x;
    suY = y;
    suVX = vx;
    suVY = vy;
    suMasa = m;
    suRadio = r;
}

void Cuerpo::paso(double dt, Cuerpo * cuerpo2) {
    double k1x = suVX;
    double k1y = suVY;
    double k1vx, k1vy;
    cuerpo2->calcular_aceleracion(suX, suY, k1vx, k1vy, suMasa);
    double k2x = suVX + (dt * k1vx / 2);
    double k2y = suVY + (dt * k1vy / 2);
    double k2vx, k2vy;
    cuerpo2->calcular_aceleracion(suX + (dt * k1x / 2), suY + (dt * k1y /2), k2vx, k2vy, suMasa);
    double k3x = suVX + (dt * k2vx / 2);
    double k3y = suVY + (dt * k2vy / 2);
    double k3vx, k3vy;
    cuerpo2->calcular_aceleracion(suX + (dt * k2x / 2), suY + (dt * k2y /2), k3vx, k3vy, suMasa);
    double k4x = suVX + dt * k3vx;
    double k4y = suVY + dt * k3vy;
    double k4vx, k4vy;
    cuerpo2->calcular_aceleracion(suX + dt * k3x, suY + dt * k3y, k4vx, k4vy, suMasa);
    suX = suX + dt * (k1x+2*k2x+2*k3x+k4x) / 6;
    suY = suY + dt * (k1y+2*k2y+2*k3y+k4y) / 6;
    suVX = suVX + dt * (k1vx+2*k2vx+2*k3vx+k4vx) / 6;
    suVY = suVY + dt * (k1vy+2*k2vy+2*k3vy+k4vy) / 6;
}

void Cuerpo::calcular_aceleracion(double x, double y, double & ax, double &ay, Cuerpo  * cuerpo2) {
    ax = -((G * (suMasa + masa2) ) * x) / ( (pow(x, 2) + pow(y, 2)) * sqrt(pow(x, 2) + pow(y, 2)));
    ay = -((G * (suMasa + masa2) ) * y) / ( (pow(x, 2) + pow(y, 2)) * sqrt(pow(x, 2) + pow(y, 2)));
}

int main() {
    /*Cuerpo sol("Sol", 0, 0, 0, 0, SOL_MASA, SOL_RADIO);
    Cuerpo tierra("Tierra", TIERRA_X0, TIERRA_Y0, TIERRA_VX0, TIERRA_VY0, TIERRA_MASA, TIERRA_RADIO);

    FILE *datos;
    datos = fopen("intento.dump","w");

    fprintf(datos,"NUM_CUERPOS\n%d\n\nNUM_PASOS\n%d\n\nNOMBRES\nSol\nTierra\n\n", NUM_CUERPOS, NUM_PASOS);
    fprintf(datos, "MASAS\n%e\n%e\n\n", SOL_MASA, TIERRA_MASA);
    fprintf(datos, "RADIOS\n%e\n%e\n\nTRAYECTORIAS\n", SOL_RADIO, TIERRA_RADIO);

    for(int nstep = 0; nstep <= 365; nstep++) {
        fprintf(datos, "%d %f %f %f %f\n", nstep, sol.x(), sol.y(), tierra.x(), tierra.y());
        tierra.paso(DT, &sol);
    }

    fclose(datos);*/

    Cuerpo pluton("Pluton", 2035706.0, 0, 230, 1, PLUTON_MASA, PLUTON_RADIO);
    Cuerpo caronte("Caronte", -17535294.0 , CARONTE_Y0, CARONTE_VX0, -230.00, CARONTE_MASA, CARONTE_RADIO);
    //Cuerpo pluton("Pluton", 0, 0, 0, 0, SOL_MASA, SOL_RADIO);
    //Cuerpo caronte("Caronte", TIERRA_X0, 0, 0, TIERRA_VY0, TIERRA_MASA, TIERRA_RADIO);
    double vxCentroDeMasa = (pluton.obtenerMasa() * pluton.vx() + caronte.obtenerMasa() * caronte.vx()) / (pluton.obtenerMasa() + caronte.obtenerMasa());
    double vyCentroDeMasa = (pluton.obtenerMasa() * pluton.vy() + caronte.obtenerMasa() * caronte.vy()) / (pluton.obtenerMasa() + caronte.obtenerMasa());
    FILE *datos;
    datos = fopen("intento.dump","w");

    fprintf(datos,"NUM_CUERPOS\n%d\n\nNUM_PASOS\n%d\n\nNOMBRES\nPluton\nCaronte\n\n", NUM_CUERPOS, NUM_PASOS);
    fprintf(datos, "MASAS\n%e\n%e\n\n", PLUTON_MASA, CARONTE_MASA);
    fprintf(datos, "RADIOS\n%e\n%e\n\nTRAYECTORIAS\n", PLUTON_RADIO, CARONTE_RADIO);

    for(int nstep = 0; nstep <= 365; nstep++) {
        caronte.asignarPosX(caronte.x()-vxCentroDeMasa * nstep * DT);
        caronte.asignarPosY(caronte.y()-vyCentroDeMasa * nstep * DT);
        pluton.asignarPosX(pluton.x()-vxCentroDeMasa * nstep * DT);
        pluton.asignarPosY(pluton.y()-vyCentroDeMasa * nstep * DT);
        fprintf(datos, "%d %f %f %f %f\n", nstep, pluton.x(), pluton.y(), caronte.x(), caronte.y());
        Cuerpo depaso("Caronte", caronte.x(),caronte.y(), caronte.vx(), caronte.vy(), caronte.obtenerMasa(), caronte.obtenerRadio());
        caronte.paso(DT, &pluton);
        pluton.paso(DT, &depaso);
    }

    fclose(datos);

    return 0;
}
