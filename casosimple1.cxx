#include <iostream>
#include <cmath>
#include <string>

#define G 1
#define SOL_MASA 1.989e30
#define SOL_RADIO 695700000.0
//#define TIERRA_X0 -147095000000.0
#define TIERRA_Y0 0.0
#define TIERRA_VX0 0.0
//#define TIERRA_VY0 -30300.0
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
#define DT 0.03
#define NUM_CUERPOS 2
#define NUM_PASOS 250*3

using namespace std;

class Resolver_r {
    public:
        Resolver_r (double x, double y, double vx, double vy, double m);
        ~Resolver_r () {}

        void asignarPosX (double x) { suX = x; }
        void asignarPosY (double y) { suY = y; }
        void asignarVelX (double x) { suVX = x; }
        void asignarVelY (double y) { suVY = y; }
        void asignarMasa (double m) { suMasa = m; }

        double x () { return suX; }
        double y () { return suY; }
        double vx () { return suVX; }
        double vy () { return suVY; }
        double masa () { return suMasa; }

        void aproximar(double x, double y, double & ax, double & ay);
        void paso(double dt);

    private:
        double suX;
        double suY;
        double suVX;
        double suVY;
        double suMasa;
};

Resolver_r::Resolver_r (double x, double y, double vx, double vy, double m) {
    suX = x;
    suY = y;
    suVX = vx;
    suVY = vy;
    suMasa = m;
}

void Resolver_r::paso(double dt) {
    double k1x = suVX;
    double k1y = suVY;
    double k1vx, k1vy;
    aproximar(suX, suY, k1vx, k1vy);
    double k2x = suVX + (dt * k1vx / 2);
    double k2y = suVY + (dt * k1vy / 2);
    double k2vx, k2vy;
    aproximar(suX + (dt * k1x / 2), suY + (dt * k1y /2), k2vx, k2vy);
    double k3x = suVX + (dt * k2vx / 2);
    double k3y = suVY + (dt * k2vy / 2);
    double k3vx, k3vy;
    aproximar(suX + (dt * k2x / 2), suY + (dt * k2y /2), k3vx, k3vy);
    double k4x = suVX + dt * k3vx;
    double k4y = suVY + dt * k3vy;
    double k4vx, k4vy;
    aproximar(suX + dt * k3x, suY + dt * k3y, k4vx, k4vy);
    suX = suX + dt * (k1x+2*k2x+2*k3x+k4x) / 6;
    suY = suY + dt * (k1y+2*k2y+2*k3y+k4y) / 6;
    suVX = suVX + dt * (k1vx+2*k2vx+2*k3vx+k4vx) / 6;
    suVY = suVY + dt * (k1vy+2*k2vy+2*k3vy+k4vy) / 6;
}

void Resolver_r::aproximar(double x, double y, double & ax, double &ay) {
    ax = -((G * suMasa) * x) / ( (pow(x, 2) + pow(y, 2)) * sqrt(pow(x, 2) + pow(y, 2)));
    ay = -((G * suMasa) * y) / ( (pow(x, 2) + pow(y, 2)) * sqrt(pow(x, 2) + pow(y, 2)));
}

int main() {

    Resolver_r solucion(2-0, 0.0-0.0, 0.01-0.05, 1.3-0.0, 1 + 1);

    double vxCentroDeMasa = (1 * (0.01) + 1 * 0.05) / solucion.masa();
    double vyCentroDeMasa = (1 * (0.5) + 1 * 0) / solucion.masa();

    FILE *datos;
    datos = fopen("casosimple1.dump","w");

    fprintf(datos,"NUM_CUERPOS\n%d\n\nNUM_PASOS\n%d\n\nNOMBRES\nPluton\nCaronte\n\n", NUM_CUERPOS, NUM_PASOS);
    fprintf(datos, "MASAS\n%e\n%e\n\n", 1, 1);
    fprintf(datos, "RADIOS\n%e\n%e\n\nTRAYECTORIAS\n", 0.5, 0.5);

    double carontex, carontey, plutonx, plutony;

    for(int nstep = 0; nstep <= NUM_PASOS; nstep++) {
        plutonx = vxCentroDeMasa * DT * nstep - solucion.x() * 1 / solucion.masa();
        plutony = vyCentroDeMasa * DT * nstep - solucion.y() * 1 / solucion.masa();
        carontex = vxCentroDeMasa * DT * nstep + solucion.x() * 1 / solucion.masa();
        carontey = vyCentroDeMasa * DT * nstep + solucion.y() * 1 / solucion.masa();
        fprintf(datos, "%d %f %f %f %f\n", nstep, plutonx, plutony, carontex, carontey);
        solucion.paso(DT);
        //solucion.asignarPosY(solucion.y()+(carontey-plutony));
    }

    fclose(datos);

    return 0;
}
