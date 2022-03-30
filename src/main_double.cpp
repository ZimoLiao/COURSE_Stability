#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <complex>
#include <math.h>
#include <omp.h>
#include <time.h>

using namespace std;

constexpr double Re = 2000, alpha = 0.179;
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
    // omp_set_num_threads(4);

    clock_t start, end;
    start = clock();

    // spatial discretization
    int ny = 2000;
    double ly = 60, dy = ly / ny;

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

    // shooting grid initialization
    int ncr = 200, nci = 200;
    double crmin = 0.0, crmax = 1.0, cimin = -0.8, cimax = 0.1;
    double *cr, *ci;
    cr = new double[ncr];
    for (int r = 0; r != ncr; r++)
        cr[r] = crmin + r * (crmax - crmin) / (ncr);
    ci = new double[nci];
    for (int i = 0; i != nci; i++)
        ci[i] = cimin + i * (cimax - cimin) / (nci - 1);

    double **Dr, **Di;
    Dr = new double *[ncr];
    Di = new double *[ncr];
    for (int r = 0; r != ncr; r++)
    {
        Dr[r] = new double[nci];
        Di[r] = new double[nci];
    }

    // solve
    {
        complex<double> p, Y[6], Yn[6],
            k1[6], k2[6], k3[6], k4[6];

        //#pragma omp parallel for
        for (int r = 0; r < ncr; r++)
        {
            for (int i = 0; i < nci; i++)
            {
                c = complex<double>(cr[r], ci[i]);
                p = sqrt(alpha * alpha + iunit * alpha * Re * (unit - c));

                Y[0] = 1;
                Y[1] = -(alpha + p);
                Y[2] = alpha * alpha + alpha * p + p * p;
                Y[3] = alpha * p;
                Y[4] = Y[3] * Y[1];
                Y[5] = Y[3] * Y[3];

                complex<double> coef = (alpha - p) * exp(-(alpha + p) * ly);
                Y[0] *= coef;
                Y[1] *= coef;
                Y[2] *= coef;
                Y[3] *= coef;
                Y[4] *= coef;
                Y[5] *= coef;

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
                }

                Dr[r][i] = Y[0].real();
                Di[r][i] = Y[0].imag();
            }
        }
    }

    end = clock();
    cout << "time = " << double(end - start) / CLOCKS_PER_SEC << "s" << endl;

    // write D
    {
        ofstream fout;
        string fname = "output/discriminate." + to_string(int(Re)) + ".dat";
        fout.open(fname);

        fout << "variables = cr ci Dr Di D lnD\n"
             << "zone        I=" + to_string(ncr) + " J=" + to_string(nci) + " F=point\n";
        for (int i = 0; i != nci; i++)
        {
            for (int r = 0; r != ncr; r++)
            {
                double D = sqrt(Dr[r][i] * Dr[r][i] + Di[r][i] * Di[r][i]);
                fout << setw(8) << cr[r] << '\t'
                     << setw(8) << ci[i] << '\t'
                     << setw(12) << Dr[r][i] << '\t'
                     << setw(12) << Di[r][i] << '\t'
                     << setw(12) << D << '\t'
                     << setw(12) << log(D + 3e-308) << '\n';
            }
        }

        fout.close();
    }

    // deallocation
    delete[] y;
    delete[] u;
    delete[] ddu;
    delete[] cr;
    delete[] ci;
    for (int r = 0; r != ncr; r++)
    {
        delete[] Dr[r];
        delete[] Di[r];
    }
    delete[] Dr;
    delete[] Di;
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