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
#define DT 1511.0
#define NUM_CUERPOS 2
#define NUM_PASOS 365

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


    //r1 pluton, r2 caronte
    double PLUTON_VY0 = 20.00;
    Resolver_r solucion(-17535294.0-2035706.0, 0.00-0.00, 0.00-0.00, -210.00-PLUTON_VY0, CARONTE_MASA + PLUTON_MASA);

    double vxCentroDeMasa = 0;
    double vyCentroDeMasa = 0*(PLUTON_MASA * (PLUTON_VY0) + CARONTE_MASA * (-210.00)) / solucion.masa();

    FILE *datos;
    datos = fopen("intento.dump","w");

    fprintf(datos,"NUM_CUERPOS\n%d\n\nNUM_PASOS\n%d\n\nNOMBRES\nPluton\nCaronte\n\n", NUM_CUERPOS, NUM_PASOS);
    fprintf(datos, "MASAS\n%e\n%e\n\n", PLUTON_MASA, CARONTE_MASA);
    fprintf(datos, "RADIOS\n%e\n%e\n\nTRAYECTORIAS\n", PLUTON_RADIO, CARONTE_RADIO);

    double carontex, carontey, plutonx, plutony;

    for(int nstep = 0; nstep <= 700; nstep++) {
        plutonx = -solucion.x() * CARONTE_MASA / solucion.masa();
        plutony = vyCentroDeMasa * DT * nstep - solucion.y() * CARONTE_MASA / solucion.masa();
        carontex = solucion.x() * PLUTON_MASA / solucion.masa();
        carontey = vyCentroDeMasa * DT * nstep + solucion.y() * PLUTON_MASA / solucion.masa();
        fprintf(datos, "%d %f %f %f %f\n", nstep, plutonx, plutony, carontex, carontey);
        solucion.paso(DT);
    }

    fclose(datos);

    /*
    r1 sol, r2 tierra
    double SOL_X0 = TIERRA_MASA * TIERRA_X0 / SOL_MASA;
    Resolver_r solucion(TIERRA_X0-SOL_X0, 0.00-0.00, 0.00-0.00, TIERRA_VY0-10.00, TIERRA_MASA + SOL_MASA);

    double vxCentroDeMasa = 0;
    double vyCentroDeMasa = (SOL_MASA * (10.00) + TIERRA_MASA * TIERRA_VY0) / solucion.masa();

    FILE *datos;
    datos = fopen("intento.dump","w");

    fprintf(datos,"NUM_CUERPOS\n%d\n\nNUM_PASOS\n%d\n\nNOMBRES\nPluton\nCaronte\n\n", NUM_CUERPOS, NUM_PASOS);
    fprintf(datos, "MASAS\n%e\n%e\n\n", PLUTON_MASA, CARONTE_MASA);
    fprintf(datos, "RADIOS\n%e\n%e\n\nTRAYECTORIAS\n", PLUTON_RADIO, CARONTE_RADIO);

    double carontex, carontey, plutonx, plutony;

    for(int nstep = 0; nstep <= 365; nstep++) {
        plutonx = -solucion.x() * TIERRA_MASA / solucion.masa();
        plutony = vyCentroDeMasa * DT * nstep - solucion.y() * TIERRA_MASA / solucion.masa();
        carontex = solucion.x() * SOL_MASA / solucion.masa();
        carontey = vyCentroDeMasa * DT * nstep + solucion.y() * SOL_MASA / solucion.masa();
        fprintf(datos, "%d %f %f %f %f\n", nstep, plutonx, plutony, carontex, carontey);
        solucion.paso(DT);
        //solucion.asignarPosY(solucion.y()+(carontey-plutony));
    }

    fclose(datos);*/

    return 0;
}
