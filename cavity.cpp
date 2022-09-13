// coding with UTF-8

// ２次元キャビティー流れのシミュレーション

#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <stdio.h>
#include <sstream>
#include <omp.h>
#include <iomanip>
using namespace std;
#define vvd vector<vector<double>>

inline void bndset(int nx, int ny, double u_wall,
                   vvd &uu, vvd &dd, vvd &pp);
inline void velocity(int nx, int ny, double dt, double dxinv, double dyinv, double dxdxinv, double dydyinv, double reinv,
                     vvd &uu, vvd &vv, vvd &pp, vvd &uup, vvd &vvp);
inline void poisson(int nx, int ny, double dtinv, double dxinv, double dyinv, double dxdxinv, double dydyinv,
                    vvd &uup, vvd &vvp, vvd &pphi);
inline void update(int nx, int ny, double dt, double dxinv, double dyinv,
                   vvd &uu, vvd &vv, vvd &uup, vvd &vvp, vvd &pp, vvd &pphi);
void output(int nx, int ny, double dx, double dy, int n,
            vvd &uu, vvd &vv, vvd &pp);

int main(int argc, char const *argv[])
{
    const int step = 100000;
    int dataout = 1000;
    int sysout = 1000;
#ifdef _OPENMP
    cout << "OpenMP threads = " << omp_get_max_threads() << endl;
#else
    cout << "No OpenMP" << endl;
#endif


    // System set ======================================== //
    const double dt = 0.0001;  // delta t
    const int nx = 256;        // x格子点数 （圧力の定義点数）
    const int ny = 256;        // y格子点数
    const double xl = 1.0;     // x代表長さ
    const double yl = 1.0;     // y代表長さ
    const double re = 100.0;   // レイノルズ数
    const double U_wall = 1.0; // 上壁流速
    // System set ======================================== //

    double dx = 1.0 / (nx - 2);
    double dy = 1.0 / (ny - 2);
    double dxinv = 1.0 / dx;
    double dyinv = 1.0 / dy;
    double dxdxinv = 1.0 / dx / dx;
    double dydyinv = 1.0 / dy / dy;
    double dtinv = 1.0 / dt;
    double reinv = 1.0 / re;

    // u[x index][y index]
    // size : u[0:nx][0:ny-1]
    vvd u(nx + 1, vector<double>(ny, 0.0));
    vvd up(nx + 1, vector<double>(ny, 0.0));

    // size : v[0:nx-1][0:ny]
    vvd v(nx, vector<double>(ny + 1, 0.0));
    vvd vp(nx, vector<double>(ny + 1, 0.0));

    // size : p[0:nx-1][0:ny-1]
    vvd p(nx, vector<double>(ny, 0.0));
    vvd phi(nx, vector<double>(ny, 0.0));

    cout << "cavity flow simulation" << endl;

    bndset(nx, ny, U_wall, u, v, p);
    bndset(nx, ny, U_wall, up, vp, phi);
    output(nx, ny, dx, dy, 0, u, v, p);

    for (int s = 1; s <= step; s++)
    {
        bndset(nx, ny, U_wall, u, v, p);                                              // u,v,p の境界条件設定
        bndset(nx, ny, U_wall, up, vp, phi);                                          // up,vp,phi の境界条件設定
        velocity(nx, ny, dt, dxinv, dyinv, dxdxinv, dydyinv, reinv, u, v, p, up, vp); // 速度の予測
        poisson(nx, ny, dtinv, dxinv, dyinv, dxdxinv, dydyinv, up, vp, phi);          // 圧力方程式を解く
        update(nx, ny, dt, dxinv, dyinv, u, v, up, vp, p, phi);                       // 更新
        if (s % dataout == 0)
            output(nx, ny, dx, dy, s, u, v, p); // output data
        if (s % sysout == 0)
            cout << s << "\ttimes done ..." << endl; // system out
    }

    return 0;
}

// boundary set
inline void bndset(int nx, int ny, double u_wall, vvd &uu, vvd &vv, vvd &pp)
{

    int i, j;
    // uu ==================================================
    // u i index
    //  0    1     2  ...... nx-2  nx-1  nx
    // buf  wall  ful ...... ful   wall  buf

    //   u j index
    //   0       buf
    // (0.5)     wall  <= moving
    //   1       ful
    //   :        :
    //   :        :
    //  ny-2     ful
    // (ny-1.5)  wall
    //  ny-1     buf

#ifdef _OPENMP
#pragma omp parallel for
#endif
    // u left + right
    for (j = 0; j < ny; j++)
    {
        // left
        uu[0][j] = uu[2][j];
        uu[1][j] = 0.0;

        // right
        uu[nx][j] = uu[nx - 2][j];
        uu[nx - 1][j] = 0.0;
    }

#ifdef _OPENMP
#pragma omp parallel for
#endif
    // u top + bottom
    for (i = 0; i < nx + 1; i++)
    {
        // top
        // (u[i][0] + u[i][1]) / 2 = u_wall
        uu[i][0] = (double)2.0 * u_wall - uu[i][1];

        // bottom
        uu[i][ny - 1] = -uu[i][ny - 2];
    }

    // vv ==================================================
    //  v i index
    //  0  (0.5)   1  ...... nx-2  (nx-1.5)  nx-1
    // buf  wall  ful ...... ful     wall     buf

    //   v j index
    //   0       buf
    //   1       wall   (<= moving)
    //   2       ful
    //   :        :
    //   :        :
    //  ny-2     ful
    //  ny-1     wall
    //   ny      buf

#ifdef _OPENMP
#pragma omp parallel for
#endif
    // v left + right
    for (j = 0; j < ny + 1; j++)
    {
        // left
        vv[0][j] = -vv[1][j];

        // right
        vv[nx - 1][j] = -vv[nx - 2][j];
    }

#ifdef _OPENMP
#pragma omp parallel for
#endif
    // v top + bottom
    for (i = 0; i < nx; i++)
    {
        // top
        vv[i][0] = vv[i][2];
        vv[i][1] = 0.0;

        // bottom
        vv[i][ny] = vv[i][ny - 2];
        vv[i][ny - 1] = 0.0;
    }

    // pp ==================================================
    //  p i index
    //  0  (0.5)   1  ...... nx-2  (nx-1.5)  nx-1
    // buf  wall  ful ...... ful     wall     buf

    //   p j index
    //   0       buf
    // (0.5)     wall  <= moving
    //   1       ful
    //   :        :
    //   :        :
    //  ny-2     ful
    // (ny-1.5)  wall
    //  ny-1     buf

#ifdef _OPENMP
#pragma omp parallel for
#endif
    // left + right
    for (j = 0; j < ny; j++)
    {
        // left
        pp[0][j] = pp[1][j];
        // right
        pp[nx - 1][j] = pp[nx - 2][j];
    }

#ifdef _OPENMP
#pragma omp parallel for
#endif
    // top + bottom
    for (i = 0; i < nx; i++)
    {
        // top
        pp[i][0] = pp[i][1];
        // bottom
        pp[i][ny - 1] = pp[i][ny - 2];
    }
}

// predict velocity
inline void velocity(int nx, int ny, double dt, double dxinv, double dyinv, double dxdxinv, double dydyinv, double reinv, vvd &uu, vvd &vv, vvd &pp, vvd &uup, vvd &vvp)
{
    double umid, vmid, udif, vdif, uad, vad;
    int i, j;
#ifdef _OPENMP
#pragma omp parallel for private(i, j, vmid, uad, udif)
#endif
    for (i = 2; i < nx - 1; i++)
    {
        for (j = 1; j < ny - 1; j++)
        {
            vmid = (vv[i - 1][j] + vv[i - 1][j + 1] + vv[i][j] + vv[i][j + 1]) * 0.25;
            uad = uu[i][j] * ((uu[i + 1][j] - uu[i - 1][j]) * 0.5 * dxinv) +
                  vmid * ((uu[i][j + 1] - uu[i][j - 1]) * 0.5 * dyinv);
            udif = (uu[i + 1][j] - 2.0 * uu[i][j] + uu[i - 1][j]) * dxdxinv +
                   (uu[i][j + 1] - 2.0 * uu[i][j] + uu[i][j - 1]) * dydyinv;
            uup[i][j] = uu[i][j] + dt * (-uad - ((pp[i][j] - pp[i - 1][j]) * dxinv) + reinv * udif);
        }
    }
#ifdef _OPENMP
#pragma omp parallel for private(i, j, umid, vad, vdif)
#endif
    for (i = 1; i < nx - 1; i++)
    {
        for (j = 2; j < ny - 1; j++)
        {
            umid = (uu[i][j - 1] + uu[i][j] + uu[i + 1][j - 1] + uu[i + 1][j]) * 0.25;
            vad = vv[i][j] * ((vv[i][j + 1] - vv[i][j - 1]) * 0.5 * dyinv) +
                  umid * ((vv[i + 1][j] - vv[i - 1][j]) * 0.5 * dxinv);
            vdif = (vv[i + 1][j] - 2.0 * vv[i][j] + vv[i - 1][j]) * dxdxinv +
                   (vv[i][j + 1] - 2.0 * vv[i][j] + vv[i][j - 1]) * dydyinv;
            vvp[i][j] = vv[i][j] + dt * (-vad - ((pp[i][j] - pp[i][j - 1]) * dyinv) + reinv * vdif);
        }
    }
}

// calculate poisson by SOR
inline void poisson(int nx, int ny, double dtinv, double dxinv, double dyinv, double dxdxinv, double dydyinv, vvd &uup, vvd &vvp, vvd &pphi)
{
    double eps = 0.00001;
    double error, divup, resid, dphi;
    double alpha = 1.6;
    long MAX = nx * ny * 2;
    long i, j;

#ifdef _OPENMP
#pragma omp parallel for private(i, j)
#endif
    for (i = 0; i < nx; i++)
    {
        for (j = 0; j < ny; j++)
        {
            pphi[i][j] = 0.0;
        }
    }

    for (long k = 0; k < MAX; k++)
    {
        error = 0.0;
#ifdef _OPENMP
#pragma omp parallel for private(i, j, divup, resid, dphi) reduction(+ \
                                                                     : error)
#endif
        for (i = 1; i < nx - 1; i++)
        {
            for (j = 1; j < ny - 1; j++)
            {
                divup = ((uup[i + 1][j] - uup[i][j]) * dxinv) +
                        ((vvp[i][j + 1] - vvp[i][j]) * dyinv);
                resid = ((pphi[i + 1][j] - 2.0 * pphi[i][j] + pphi[i - 1][j]) * dxdxinv) +
                        ((pphi[i][j + 1] - 2.0 * pphi[i][j] + pphi[i][j - 1]) * dydyinv) -
                        dtinv * divup;
                dphi = alpha * resid / (2.0 * (dxdxinv + dydyinv));
                error += abs(dphi);
                pphi[i][j] = pphi[i][j] + dphi;
            }
        }
        if ((error / (nx * ny)) < eps)
            break;
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (j = 0; j < ny; j++)
        {
            pphi[0][j] = pphi[1][j];
            pphi[nx - 1][j] = pphi[nx - 2][j];
        }
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (i = 0; i < nx; i++)
        {
            pphi[i][0] = pphi[i][1];
            pphi[i][ny - 1] = pphi[i][ny - 2];
        }
    }
}

// update new u,v,p
inline void update(int nx, int ny, double dt, double dxinv, double dyinv, vvd &uu, vvd &vv, vvd &uup, vvd &vvp, vvd &pp, vvd &pphi)
{
    int i, j;
#ifdef _OPENMP
#pragma omp parallel for private(j)
#endif
    for (i = 2; i < nx - 1; i++)
    {
        for (j = 1; j < ny - 1; j++)
        {
            uu[i][j] = uup[i][j] - dt * (pphi[i + 1][j] - pphi[i - 1][j]) * 0.5 * dxinv;
        }
    }
#ifdef _OPENMP
#pragma omp parallel for private(j)
#endif
    for (i = 1; i < nx - 1; i++)
    {
        for (j = 2; j < ny - 1; j++)
        {
            vv[i][j] = vvp[i][j] - dt * (pphi[i][j + 1] - pphi[i][j - 1]) * 0.5 * dyinv;
        }
    }
#ifdef _OPENMP
#pragma omp parallel for private(j)
#endif
    for (i = 1; i < nx - 1; i++)
    {
        for (j = 1; j < ny - 1; j++)
        {
            pp[i][j] = pp[i][j] + pphi[i][j];
        }
    }
}

// output data
void output(int nx, int ny, double dx, double dy, int n, vvd &uu, vvd &vv, vvd &pp)
{
    ostringstream zfill;
    zfill << setfill('0') << setw(6) << n;
    string fname = "data/xyuvp" + zfill.str() + ".tsv";
    ofstream ofs(fname);
    if (!ofs)
    {
        cerr << " ERROR : cannot open " << fname << endl;
        exit(-1);
    }

    for (int i = 1; i < nx - 1; i++)
    {
        for (int j = 1; j < ny - 1; j++)
        {
            ofs << (double)i * dx << '\t';
            ofs << (double)j * dy << '\t';
            ofs << 0.5 * (uu[i][j] + uu[i + 1][j]) << '\t';
            ofs << 0.5 * (vv[i][j] + vv[i][j + 1]) << '\t';
            ofs << pp[i][j] << '\n';
        }
        ofs << '\n';
    }
    ofs.close();
}
