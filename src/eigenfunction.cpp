#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <complex>
#include <math.h>
#include <omp.h>
#include <time.h>

using namespace std;

constexpr double Re = 580, alpha = 0.179;
constexpr complex<double> iunit = {0.0, 1.0};
double unit = 1.0;

// r.h.s. of blasius equation
void f_blasius(const double var[3],
               double f[3]);

// r.h.s. of compound matrix for O-S equation
void f_compound(const complex<double> c, const double u, const double ddu,
                const complex<double> var[6],
                complex<double> f[6]);

int main()
{
    // spatial discretization
    int ny = 5000;
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
        double var[3] = {0.0, 0.0, 0.332057336}, varn[3], // TODO: Blasius solution
            k1[3], k2[3], k3[3], k4[3];

        for (int i = 0; i < ny - 1; i++)
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
            ddu[i + 1] = -var[0] * var[2] / 2.0;
        }
    }

    // TODO: modify
    // for (int i = 0; i != ny; i++)
    // {
    //     u[i] /= u[ny - 1];
    // }

    // write blasius profile
    {
        ofstream fout;
        fout.open("output/blasius.dat");

        fout << "variables = y u ddu\n";
        for (int i = 0; i != ny; i++)
        {
            fout << setw(8) << y[i] << '\t' << setw(12) << u[i] << '\t' << setw(12) << ddu[i] << '\n';
        }

        fout.close();
    }

    /* compound matrix - shooting method */
    // parameters
    // double Re = 580, alpha = 0.179;
    complex<double> c;
    double cr = 0.36412269, ci = 0.00795979;
    double Dr, Di;
    double Yr[10000][6], Yi[10000][6];

    // solve
    {
        complex<double> p, Y[6], Yn[6],
            k1[6], k2[6], k3[6], k4[6];

        c = complex<double>(cr, ci);
        p = sqrt(alpha * alpha + iunit * alpha * Re * (unit - c));

        Y[0] = 1;
        Y[1] = -(alpha + p);
        Y[2] = alpha * alpha + alpha * p + p * p;
        Y[3] = alpha * p;
        Y[4] = Y[3] * Y[1];
        Y[5] = Y[3] * Y[3];

        complex<double> coef = (alpha - p) * exp(-(alpha + p) * ly);
        for (int v = 0; v != 6; v++)
        {
            Y[v] *= coef;
            Yr[ny][v] = Y[v].real();
            Yi[ny][v] = Y[v].imag();
        }

        // cout << r << '\t' << i << '\t' << abs(coef) << '\n';

        for (int step = ny - 1; step >= 0; step--)
        {
            f_compound(c, u[step], ddu[step], Y, k1);
            // if (step != 0)
            // {
            //     for (int v = 0; v != 6; v++)
            //         Yn[v] = Y[v] - k1[v] * dy / 2.0;
            //     f_compound(c, (u[step] + u[step - 1]) / 2.0, (ddu[step] + ddu[step - 1]) / 2.0, Yn, k2);

            //     for (int v = 0; v != 6; v++)
            //         Yn[v] = Y[v] - k2[v] * dy / 2.0;
            //     f_compound(c, (u[step] + u[step - 1]) / 2.0, (ddu[step] + ddu[step - 1]) / 2.0, Yn, k3);

            //     for (int v = 0; v != 6; v++)
            //         Yn[v] = Y[v] - k4[v] * dy;
            //     f_compound(c, u[step - 1], ddu[step - 1], Yn, k4);

            //     for (int v = 0; v != 6; v++)
            //         Y[v] -= (k1[v] + 2.0 * k2[v] + 2.0 * k3[v] + k4[v]) * dy / 6.0;
            // }
            // else
            if (step != 0) // RK-3
            {
                for (int v = 0; v != 6; v++)
                    Yn[v] = Y[v] - k1[v] * dy / (unit * 2.0);
                f_compound(c, (u[step] + u[step - 1]) / 2.0, (ddu[step] + ddu[step - 1]) / 2.0, Yn, k2);

                for (int v = 0; v != 6; v++)
                    Yn[v] = Y[v] - (unit * 2.0) * k2[v] * dy + k1[v] * dy;
                f_compound(c, u[step - 1], ddu[step - 1], Yn, k3);

                for (int v = 0; v != 6; v++)
                    Y[v] -= (k1[v] + (unit * 4.0) * k2[v] + k3[v]) * dy / (unit * 6.0);
            }
            else
            {
                for (int v = 0; v != 6; v++)
                    Y[v] -= k1[v] * dy;
            }

            for (int v = 0; v != 6; v++)
            {
                Yr[step][v] = Y[v].real();
                Yi[step][v] = Y[v].imag();
            }
        }

        Dr = Y[0].real();
        Di = Y[0].imag();
    }

    cout << Dr << '\t' << Di << '\n';

    // write Y
    ofstream fout;
    fout.open("output/Y.dat");
    fout << "variables = eta y1r y1i y2r y2i y3r y3i y4r y4i y5r y5i y6r y6i\n";
    for (int i = 0; i != ny; i++)
    {
        fout << i * dy << '\t';
        for (int v = 0; v != 6; v++)
        {
            fout << Yr[i][v] << '\t' << Yi[i][v] << '\t';
        }
        fout << '\n';
    }
    fout.close();

    // calculate eigenfunction (using eqn.3)
    double phir[10000] = {0.0}, phii[10000] = {0.0}, phimax = 0.0;
    {
        complex<double> phi[3] = {0.0, 0.0, (1.0 - iunit) * sqrt(Re)}, dphi[3];
        for (int step = 0; step != ny; step++)
        {
            dphi[0] = phi[1];
            dphi[1] = phi[2];
            dphi[2] = (complex<double>(Yr[step][2], Yi[step][2]) * phi[2] - complex<double>(Yr[step][5], Yi[step][5]) * phi[0]) / complex<double>(Yr[step][1], Yi[step][1]);

            for (int v = 0; v != 3; v++)
            {
                phi[v] += dphi[v] * dy;
            }

            phir[step + 1] = phi[1].real();
            phii[step + 1] = phi[1].imag();
            if (phimax < sqrt(phir[step + 1] * phir[step + 1] + phii[step + 1] * phii[step + 1]))
            {
                phimax = sqrt(phir[step + 1] * phir[step + 1] + phii[step + 1] * phii[step + 1]);
            }
        }
    }

    // write eigenfunction
    fout.open("output/eigenfunction.dat");
    fout << "variables = eta phir phii phi\n";
    for (int i = 0; i <= ny; i+=8)
    {
        fout << i * dy << '\t' << phir[i] / phimax << '\t' << phii[i] / phimax << '\t' << sqrt(phir[i] * phir[i] + phii[i] * phii[i]) / phimax << '\n';
    }
    fout.close();

    // deallocation
    delete[] y;
    delete[] u;
    delete[] ddu;
}

void f_blasius(const double var[3],
               double f[3])
{
    f[0] = var[1];
    f[1] = var[2];
    f[2] = -var[0] * var[2] / 2.0;
}

void f_compound(const complex<double> c, const double u, const double ddu,
                const complex<double> var[6],
                complex<double> f[6])
{
    complex<double> a2, a4;
    a2 = 2.0 * alpha * alpha + iunit * alpha * Re * (u - c);
    a4 = -pow(alpha, 4.0) - iunit * alpha * Re * (alpha * alpha * (u - c) + ddu);

    f[0] = var[1];
    f[1] = var[2] + var[3];
    f[2] = a2 * var[1] + var[4];
    f[3] = var[4];
    f[4] = -a4 * var[0] + a2 * var[3] + var[5];
    f[5] = -a4 * var[1];
}