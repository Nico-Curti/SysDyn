#define SAVEDAT
#include <iostream>
#include <random>
#include <algorithm>
#include <cstring>
#ifdef VIEWER
// g++ diffusion2D.cpp -std=c++11 `pkg-config opencv --cflags --libs`
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/contrib/contrib.hpp>
#endif
#ifdef SAVEDAT
#include <fstream>
#endif
#define A 4.5f
#define B 6.75f
#define Du 2.0f
#define Dv 16.f
#define ITER 60000
#define ht .005f
#define hx 1.f
#define SEED 123

constexpr float hx2 = 1.f / hx/hx;
std::mt19937 eng(SEED); 
std::uniform_real_distribution<float> distr(0.f, 1.f);
/*
                            Brusselator - Diffusion Model

                                        A -> X
                                    B + X -> Y + C
                                   2X + Y -> 3X
                                        X -> D

We assume that A, B, C and D concetrations are constants, so we can model the
system evolution as a cinetic model reactions

1 reaction :  dX/dt = +V1    dY/dt = 0
2 reaction :  dX/dt = -V2    dY/dt = +V2
3 reaction :  dX/dt = +V3    dY/dt = -V3
4 reaction :  dX/dt = -V4    dY/dt = 0

Considering positive and negative terms of each equation and using the Law of Mass Action
we rewrite the equations as
                        
                                dX/dt = A - BX + X^2Y - X
                                dY/dt = BX - X^2Y



*/

inline float Vx(const float &x, const float &y){ return A - (B+1)*x + y*x*x; }
inline float Vy(const float &x, const float &y){ return B*x - x*x*y; }

void diffusion2D(float **ut, float **vt, const int &dim, const float &dx)
{
    float **u = new float*[dim],
          **v = new float*[dim],
          lap_u, lap_v;
    std::generate(u, u + dim, [&dim](){return new float[dim];});
    std::generate(v, v + dim, [&dim](){return new float[dim];});
    int i, j;
    for(int t = 0; t < ITER; ++t)
    {
#ifdef DEBUG
        std::cout << "Iteration " << t << " of " << ITER << std::endl;
#endif
        for(i = 0; i < dim; ++i)
        {
            std::memcpy(u[i], ut[i], sizeof(float)*dim);
            std::memcpy(v[i], vt[i], sizeof(float)*dim);
        }
        // i == 0
        // j == 0
        lap_u = (u[0][1] + u[1][0] + u[0][dim-1] + u[dim-1][0] - 4*u[0][0]) * hx2;
        lap_v = (v[0][1] + v[1][0] + v[0][dim-1] + v[dim-1][0] - 4*v[0][0]) * hx2;
        ut[0][0] = ht * (Du * lap_u + Vx(u[0][0], v[0][0])) + u[0][0];
        vt[0][0] = ht * (Du * lap_v + Vy(u[0][0], v[0][0])) + v[0][0];
        // j == dim - 1
        lap_u = (u[0][dim-2] + u[1][dim-1] + u[0][0] + u[dim-1][0] - 4*u[0][dim-1]) * hx2;
        lap_v = (v[0][dim-2] + v[1][dim-1] + v[0][0] + v[dim-1][0] - 4*v[0][dim-1]) * hx2;
        ut[0][dim-1] = ht * (Du * lap_u + Vx(u[0][dim-1], v[0][dim-1])) + u[0][dim-1];
        vt[0][dim-1] = ht * (Dv * lap_v + Vy(u[0][dim-1], v[0][dim-1])) + v[0][dim-1];

        for(i = 1; i < dim - 1; ++i)
        {
            // first row
            lap_u = (u[0][i+1] + u[0][i-1] + u[dim-1][i] + u[1][i] - 4*u[0][i]) * hx2;
            lap_v = (v[0][i+1] + v[0][i-1] + v[dim-1][i] + v[1][i] - 4*v[0][i]) * hx2;
            ut[0][i] = ht * (Du * lap_u + Vx(u[0][i], v[0][i])) + u[0][i];
            vt[0][i] = ht * (Dv * lap_v + Vy(u[0][i], v[0][i])) + v[0][i];

            // j == 0
            lap_u = (u[i][dim-1] + u[i][1] + u[i-1][0] + u[i+1][0] - 4*u[i][0]) * hx2;
            lap_v = (v[i][dim-1] + v[i][1] + v[i-1][0] + v[i+1][0] - 4*v[i][0]) * hx2;
            ut[i][0] = ht * (Du * lap_u + Vx(u[i][0], v[i][0])) + u[i][0];
            vt[i][0] = ht * (Dv * lap_v + Vy(u[i][0], v[i][0])) + v[i][0];
            
            for(j = 1; j < dim - 1; ++j)
            {
                lap_u = (u[i][j+1] + u[i][j-1] + u[i-1][j] + u[i+1][j] - 4*u[i][j]) * hx2;
                lap_v = (v[i][j+1] + v[i][j-1] + v[i-1][j] + v[i+1][j] - 4*v[i][j]) * hx2;

                ut[i][j] = ht * (Du * lap_u + Vx(u[i][j], v[i][j])) + u[i][j];
                vt[i][j] = ht * (Dv * lap_v + Vy(u[i][j], v[i][j])) + v[i][j];
            }
            // j == dim - 1
            lap_u = (u[i][j-1] + u[i-1][j] + u[i+1][j] + u[i][0] - 4*u[i][j]) * hx2;
            lap_v = (v[i][j-1] + v[i-1][j] + v[i+1][j] + v[i][0] - 4*v[i][j]) * hx2;
            ut[i][j] = ht * (Du * lap_u + Vx(u[i][j], v[i][j])) + u[i][j];
            vt[i][j] = ht * (Dv * lap_v + Vy(u[i][j], v[i][j])) + v[i][j];

            // last row
            lap_u = (u[dim-1][i+1] + u[dim-1][i-1] + u[dim-1][i] + u[0][i] - 4*u[dim-1][i]) * hx2;
            lap_v = (v[dim-1][i+1] + v[dim-1][i-1] + v[dim-1][i] + v[0][i] - 4*v[dim-1][i]) * hx2;
            ut[dim-1][i] = ht * (Du * lap_u + Vx(u[dim-1][i], v[dim-1][i])) + u[dim-1][i];
            vt[dim-1][i] = ht * (Dv * lap_v + Vy(u[dim-1][i], v[dim-1][i])) + v[dim-1][i];
        }

        // i == dim-1
        // j == 0
        lap_u = (u[dim-1][1] + u[0][0] + u[dim-1][dim-1] + u[dim-2][0] - 4*u[dim-1][0]) * hx2;
        lap_v = (v[dim-1][1] + v[0][0] + v[dim-1][dim-1] + v[dim-2][0] - 4*v[dim-1][0]) * hx2;
        ut[dim-1][0] = ht * (Du * lap_u + Vx(u[dim-1][0], v[dim-1][0])) + u[dim-1][0];
        vt[dim-1][0] = ht * (Du * lap_v + Vy(u[dim-1][0], v[dim-1][0])) + v[dim-1][0];

        // j == dim - 1
        lap_u = (u[dim-1][dim-2] + u[dim-2][dim-1] + u[dim-1][0] + u[0][dim-1] - 4*u[dim-1][dim-1]) * hx2;
        lap_v = (v[dim-1][dim-2] + v[dim-2][dim-1] + v[dim-1][0] + v[0][dim-1] - 4*v[dim-1][dim-1]) * hx2;
        ut[dim-1][dim-1] = ht * (Du * lap_u + Vx(u[dim-1][dim-1], v[dim-1][dim-1])) + u[dim-1][dim-1];
        vt[dim-1][dim-1] = ht * (Dv * lap_v + Vy(u[dim-1][dim-1], v[dim-1][dim-1])) + v[dim-1][dim-1];
    }
#ifdef SAVEDAT
    std::ofstream os("end.dat", std::ios::out | std::ios::binary);
    for(int i = 0; i < dim; ++i)
        for(int j = 0; j < dim; ++j) 
            os.write( (const char *) &ut[i][j], sizeof( float ));
    os.close();
#endif

    for(int i = 0; i < dim; ++i)
    {
        delete[] u[i];
        delete[] v[i];
    }
    delete[] u;
    delete[] v;
    return;
}


int main(int argc, char **argv)
{
    const int dim = 120;
    const float dx = 2.f / dim;
    float **u = new float*[dim],
          **v = new float*[dim];
    std::generate(u, u + dim, [&dim](){return new float[dim];});
    std::generate(v, v + dim, [&dim](){return new float[dim];});

    // initial condition
    for(int i = 0; i < dim; ++i)
    {
        std::generate(u[i], u[i] + dim, [](){return A   + distr(eng)*.3f;});
        std::generate(v[i], v[i] + dim, [](){return B/A + distr(eng)*.3f;});
    }
#ifdef SAVEDAT
    std::ofstream os("init.dat", std::ios::out | std::ios::binary);
    for(int i = 0; i < dim; ++i)
        for(int j = 0; j < dim; ++j) 
            os.write( (const char *) &u[i][j], sizeof( float ));
    os.close();
#endif

    diffusion2D(u, v, dim, dx);
    
    return 0;
}