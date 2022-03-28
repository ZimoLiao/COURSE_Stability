#include <iostream>
#include <fstream>
#include <string>
#include <math.h>

using namespace std;

// r.h.s. of blasius equation
void f_blasius(const double var[3],
               double varn[3]);

int main()
{
    // spatial discretization
    int ny = 200;
    double ly = 10, dy = ly / ny;

    double *y = new double[ny];
    double *u = new double[ny];
    double *ddu = new double[ny];
    for (int i = 0; i != ny; i++)
    {
        y[i] = dy * i;
        u[i] = 0.0;
        ddu[i] = 0.0;
    }

    // Blasius solution (RK-4)
    {
        double var[3] = {0.0, 0.0, 0.46}, varn[3], // TODO: Blasius solution
            k1[3], k2[3], k3[3], k4[3];
        for (int i = 0; i != ny - 1; i++)
        {
            f_blasius(var, k1);

            for (int v = 0; v != 3; v++)
                varn[v] = var[v] + 0.5 * dy * k1[v];
            f_blasius(varn, k2);

            for (int v = 0; v != 3; v++)
                varn[v] = var[v] + 0.5 * dy * k2[v];
            f_blasius(varn, k3);

            for (int v = 0; v != 3; v++)
                varn[v] = var[v] + dy * k3[v];
            f_blasius(varn, k4);

            for (int v = 0; v != 3; v++)
                var[v] += dy * (k1[v] + 2.0 * k2[v] + 2.0 * k3[v] + k4[v]) / 6.0;

            u[i + 1] = var[1];
            ddu[i + 1] = -var[0] * var[2];
        }
    }

    // write blasius profile
    {
        ofstream fout;
        fout.open("blasius.dat");

        for (int i = 0; i != ny; i++)
        {
            fout << u[i] << '\t' << ddu[i] << '\n';
        }

        fout.close();
    }

    delete[] y;
    delete[] u;
}

void f_blasius(const double var[3],
               double varn[3])
{
    varn[0] = var[1];
    varn[1] = var[2];
    varn[2] = -var[0] * var[2];
}