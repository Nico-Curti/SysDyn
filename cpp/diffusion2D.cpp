#include <iostream>
#include <random>
#include <algorithm>
#include <array>
#include <memory>
#include <chrono>
#ifdef VIEWER
// g++ diffusion2D.cpp -O3 -std=c++11 `pkg-config opencv --cflags --libs`
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
//#include <opencv2/contrib/contrib.hpp>
#endif
#ifdef SAVEDAT
#include <fstream>
#endif
static constexpr float A  = 4.5f;
static constexpr float B  = 6.75f;
static constexpr float Du = 2.0f;
static constexpr float Dv = 16.f;
static constexpr int ITER = 60000;
static constexpr float ht = .005f;
static constexpr float hx = 1.f;
static constexpr int SEED = 123;

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

auto Vx = [](const float &x, const float &y){ return A - (B+1)*x + y*x*x; };
auto Vy = [](const float &x, const float &y){ return B*x - x*x*y; };

void diffusion2D(float *ut, float *vt, const int &dim, const float &dx)
{
  int i, j, dim2 = dim*dim;
  float lap_u, lap_v;
  std::unique_ptr<float[]> u(new float[dim2]),
                           v(new float[dim2]);

  for(int t = 0; t < ITER; ++t)
  {
#ifdef DEBUG
    std::cout << "Iteration " << t << " of " << ITER << std::endl;
#endif
    std::copy_n(ut, dim2, u.get());
    std::copy_n(vt, dim2, v.get());

    // i == 0
    // j == 0
    lap_u = (u[1] + u[dim] + u[dim-1] + u[dim2 - dim] - 4*u[0]) * hx2;
    lap_v = (v[1] + v[dim] + v[dim-1] + v[dim2 - dim] - 4*v[0]) * hx2;
    ut[0] = ht * (Du * lap_u + Vx(u[0], v[0])) + u[0];
    vt[0] = ht * (Du * lap_v + Vy(u[0], v[0])) + v[0];
    // j == dim - 1
    lap_u = (u[dim-2] + u[dim + dim-1] + u[0] + u[dim2-1] - 4*u[dim-1]) * hx2;
    lap_v = (v[dim-2] + v[dim + dim-1] + v[0] + v[dim2-1] - 4*v[dim-1]) * hx2;
    ut[dim-1] = ht * (Du * lap_u + Vx(u[dim-1], v[dim-1])) + u[dim-1];
    vt[dim-1] = ht * (Dv * lap_v + Vy(u[dim-1], v[dim-1])) + v[dim-1];

    for(i = 1; i < dim - 1; ++i)
    {
      // first row
      lap_u = (u[i+1] + u[i-1] + u[dim2 - dim + i] + u[dim + i] - 4*u[i]) * hx2;
      lap_v = (v[i+1] + v[i-1] + v[dim2 - dim + i] + v[dim + i] - 4*v[i]) * hx2;
      ut[i] = ht * (Du * lap_u + Vx(u[i], v[i])) + u[i];
      vt[i] = ht * (Dv * lap_v + Vy(u[i], v[i])) + v[i];

      // j == 0
      lap_u = (u[dim*i + dim-1] + u[i*dim + 1] + u[(i-1)*dim] + u[(i+1)*dim] - 4*u[dim*i]) * hx2;
      lap_v = (v[dim*i + dim-1] + v[i*dim + 1] + v[(i-1)*dim] + v[(i+1)*dim] - 4*v[dim*i]) * hx2;
      ut[dim*i] = ht * (Du * lap_u + Vx(u[dim*i], v[dim*i])) + u[dim*i];
      vt[dim*i] = ht * (Dv * lap_v + Vy(u[dim*i], v[dim*i])) + v[dim*i];

      for(j = 1; j < dim - 1; ++j)
      {
        lap_u = (u[dim*i + j+1] + u[i*dim + j-1] + u[(i-1)*dim + j] + u[(i+1)*dim + j] - 4*u[dim*i + j]) * hx2;
        lap_v = (v[dim*i + j+1] + v[i*dim + j-1] + v[(i-1)*dim + j] + v[(i+1)*dim + j] - 4*v[dim*i + j]) * hx2;

        ut[i*dim + j] = ht * (Du * lap_u + Vx(u[i*dim + j], v[i*dim + j])) + u[i*dim + j];
        vt[i*dim + j] = ht * (Dv * lap_v + Vy(u[i*dim + j], v[i*dim + j])) + v[i*dim + j];
      }
      // j == dim - 1
      lap_u = (u[i*dim + j-1] + u[(i-1)*dim + j] + u[(i+1)*dim + j] + u[i*dim] - 4*u[dim*i + j]) * hx2;
      lap_v = (v[i*dim + j-1] + v[(i-1)*dim + j] + v[(i+1)*dim + j] + v[i*dim] - 4*v[dim*i + j]) * hx2;
      ut[i*dim + j] = ht * (Du * lap_u + Vx(u[i*dim + j], v[i*dim + j])) + u[i*dim + j];
      vt[i*dim + j] = ht * (Dv * lap_v + Vy(u[i*dim + j], v[i*dim + j])) + v[i*dim + j];

      // last row
      lap_u = (u[dim2 - dim + i+1] + u[dim2 - dim + i-1] + u[dim2 - dim + i] + u[i] - 4*u[dim2 - dim + i]) * hx2;
      lap_v = (v[dim2 - dim + i+1] + v[dim2 - dim + i-1] + v[dim2 - dim + i] + v[i] - 4*v[dim2 - dim + i]) * hx2;
      ut[dim2 - dim + i] = ht * (Du * lap_u + Vx(u[dim2 - dim + i], v[dim2 - dim + i])) + u[dim2 - dim + i];
      vt[dim2 - dim + i] = ht * (Dv * lap_v + Vy(u[dim2 - dim + i], v[dim2 - dim + i])) + v[dim2 - dim + i];
    }

    // i == dim-1
    // j == 0
    lap_u = (u[dim2 - dim + 1] + u[0] + u[dim2 - 1] + u[dim2 - 2*dim] - 4*u[dim2 - dim]) * hx2;
    lap_v = (v[dim2 - dim + 1] + v[0] + v[dim2 - 1] + v[dim2 - 2*dim] - 4*v[dim2 - dim]) * hx2;
    ut[dim2 - dim] = ht * (Du * lap_u + Vx(u[dim2 - dim], v[dim2 - dim])) + u[dim2 - dim];
    vt[dim2 - dim] = ht * (Du * lap_v + Vy(u[dim2 - dim], v[dim2 - dim])) + v[dim2 - dim];

    // j == dim - 1
    lap_u = (u[dim2 - 2] + u[dim2 - dim - 1] + u[dim2 - dim] + u[dim - 1] - 4*u[dim2 - 1]) * hx2;
    lap_v = (v[dim2 - 2] + v[dim2 - dim - 1] + v[dim2 - dim] + v[dim - 1] - 4*v[dim2 - 1]) * hx2;
    ut[dim2 - 1] = ht * (Du * lap_u + Vx(u[dim2 - 1], v[dim2 - 1])) + u[dim2 - 1];
    vt[dim2 - 1] = ht * (Dv * lap_v + Vy(u[dim2 - 1], v[dim2 - 1])) + v[dim2 - 1];
  }

  return;
}


int main(int argc, char **argv)
{
  const int dim = 120;
  const float dx = 2.f / dim;
  std::array<float, dim*dim> u, v;

  // initial condition
  std::generate_n(u.begin(), dim*dim, [](){return A   + distr(eng)*.3f;});
  std::generate_n(v.begin(), dim*dim, [](){return B/A + distr(eng)*.3f;});

#ifdef SAVEDAT
  std::ofstream os("init.dat", std::ios::out | std::ios::binary);
  for(int i = 0; i < dim*dim; ++i)
    os.write( (const char *) &u[i], sizeof( float ));
  os.close();
#endif
  auto start = std::chrono::steady_clock::now();
  diffusion2D(u.data(), v.data(), dim, dx);
  auto duration = std::chrono::duration_cast<std::chrono::seconds> (std::chrono::steady_clock::now() - start);
  std::cout << "Elapsed time : " << duration.count() << " sec" << std::endl;

#ifdef SAVEDAT
  std::ofstream os("end.dat", std::ios::out | std::ios::binary);
  for(int i = 0; i < dim*dim; ++i)
    os.write( (const char *) &u[i], sizeof( float ));
  os.close();
#endif
#ifdef VIEWER
  cv::Mat U(dim, dim, CV_32FC1, u.data());
  cv::normalize(U, U, 0, 255, cv::NORM_MINMAX);
  U.convertTo( U, CV_8UC1 );
  // you need OpenCV contrib to applyColorMap !!!
  //applyColorMap(U, U, cv::COLORMAP_JET);
  cv::namedWindow( "Turing Pattern", cv::WINDOW_NORMAL );// Create a window for display
  cv::imshow( "Turing Pattern", U );
  cv::waitKey(0);
  U.release();
#endif
  return 0;
}