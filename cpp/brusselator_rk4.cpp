// g++ brusselator_rk4.cpp -std=c++14 -O3 -o brusselator_rk4

#include <memory>
#include <algorithm>
#include <cmath>
#include <type_traits>
#include <cassert>

template < class type >
using is_floating_point = typename std :: enable_if < std :: is_floating_point < type > :: value > :: type *;

template < class type, is_floating_point < type > = nullptr >
using array = std :: unique_ptr < type[] >;


/**
* @brief Brusselator kinetic
*
* @param t List of time points.
* @param x0 Initial condition of the x signal.
* @param y0 Initial condition of the y signal.
* @param A Constant of the reaction.
* @param B Constant of the reaction.
* @param x The resulting x array.
* @param y The resulting y array.
*
* @tparam type Data-type of arrays
* @tparam Length of time points.
*
*/
template < class type, int32_t N >
void Brusselator (const array < type > & t,
                  const type & x0, const type & y0, const type & A, const type & B,
                  array < type > & x, array < type > & y
                  )
{
  // determine the interval as diff
  const type dt = t[1] - t[0]; // Note: we are assuming it is constant!!

  // Set the initial condition
  x[0] = x0;
  y[0] = y0;

  // set the equation functions
  auto dx = [](const type & x, const type & y, const type & A, const type & B)
            {
              return A + x*x*y - B*x - x;
            };
  auto dy = [](const type & x, const type & y, const type & A, const type & B)
            {
              return B*x - x*x*y;
            };


  // Integrate the equations using the RK4 method
  for (int32_t i = 0; i < N - 1; ++i)
  {
    const type kx1  = dt * dx( x[i], y[i], A, B);
    const type ky1  = dt * dy( x[i], y[i], A, B);
    
    const type kx2  = dt * dx( x[i] + .5 * kx1, y[i] + .5 * ky1, A, B);
    const type ky2  = dt * dy( x[i] + .5 * kx1, y[i] + .5 * ky1, A, B);
    
    const type kx3  = dt * dx( x[i] + .5 * kx2, y[i] + .5 * ky2, A, B);
    const type ky3  = dt * dy( x[i] + .5 * kx2, y[i] + .5 * ky2, A, B);

    const type kx4  = dt * dx( x[i] + kx3, y[i] + ky3, A, B);
    const type ky4  = dt * dy( x[i] + kx3, y[i] + ky3, A, B);

    x[i + 1]  = x[i]  + type(1. / 6.) * (kx1 + type(2.) * kx2  + type(2.) * kx3 + kx4);
    y[i + 1]  = y[i]  + type(1. / 6.) * (ky1 + type(2.) * ky2  + type(2.) * ky3 + ky4);
  }
}


int32_t main (/*int32_t argc, char ** argv*/)
{
  const float x0  = 1.6f;
  const float y0  = 2.8f;

  const float A = .5f;
  const float B = 2.f;

  const int32_t iterations = 10000;
  const float dt = 1e-2f;

  array < float > time(new float[iterations]);
  std :: generate_n(time.get(), iterations, [n = 0, dt] () mutable { return dt * n++; });

  array < float > x(new float[iterations]);
  array < float > y(new float[iterations]);

  Brusselator < float, iterations >(time, x0, y0, A, B, x, y);;

  return 0;
}