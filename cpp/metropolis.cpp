// g++ metropolis.cpp -std=c++14 -O3 -o metropolis -fopenmp

#include <random>
#include <iostream>
#include <iomanip>
#include <numeric>
#include <omp.h>

int main()
{
  std :: random_device rd;
  std :: mt19937 eng{rd()};
  
  int64_t N = 1e5, N_run = 1e4;
  
  double pi = 0.;
  double start_time, run_time;
    
  std :: uniform_real_distribution < double > uniform_dist {-1., 1.};

  omp_set_num_threads(omp_get_max_threads());
  
  start_time = omp_get_wtime();
  #pragma omp parallel
  {
    #pragma omp for reduction (+ : pi)
    for (int64_t I = 0; I < N_run; ++I)
    {
      int64_t Nhints = 0;
            
      #pragma omp parallel
      {
        #pragma omp for reduction (+ : Nhints)
        for (int64_t i = 0; i < N; ++i)
        {
          const double x = uniform_dist(eng);
          const double y = uniform_dist(eng);
          const double res = x*x + y*y;
          Nhints += res < 1 ? 1 : 0;
        }
      }
      pi += 4.0 * (static_cast < double >(Nhints) / N);
    }
  }
  
  pi /= N_run;
  run_time = omp_get_wtime() - start_time;
  std :: cout << "pi with " << N << " steps for " << N_run << " runs is " << std :: setprecision(16) << pi << " in " << run_time << " sec" << std :: endl;
    
  return 0;
}   
