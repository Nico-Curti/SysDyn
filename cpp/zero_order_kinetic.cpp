// g++ zero_order_kinetic.cpp -std=c++14 -O3 -o zero_order_kinetic

#include <memory>
#include <algorithm>
#include <type_traits>

template < class type >
using is_floating_point = typename std :: enable_if < std :: is_floating_point < type > :: value > :: type *;

template < class type, is_floating_point < type > = nullptr >
using array = std :: unique_ptr < type[] >;

/**
* @brief Zero order kinetic
*
* @param x List of time points.
* @param y0 Initial condition of the reactant.
* @param alpha Constant of the reaction.
*
* @tparam type Data-type of arrays
* @tparam Length of time points.
*
* @return The resulting product array.
*
*/
template < class type, int32_t N >
array < type > zero_order (const array < type > & x, const type & y0, const type & alpha)
{
  // determine the interval as diff
  const type dx = x[1] - x[0]; // Note: we are assuming it is constant!!

  // Create an empyt buffer to store our results
  // Its length must be greater than x since we want to set
  // the initial condition!
  array < type > y = std :: make_unique < type[] >(N + 1);
  // Set the initial condition
  y[0] = y0;

  // Integrate the equation using the Euler method
  for (int32_t i = 0; i < N; ++i)
    y[i + 1] = y[i] * (type(1.) - alpha * dx);

  return y;
}


int32_t main (/*int32_t argc, char ** argv*/)
{
  const float y0 = 10.f;
  const float alpha = 0.5f;

  const int32_t iterations = 1000;
  const float dt = 1e-2f;

  array < float > x(new float[iterations]);
  std :: generate_n(x.get(), iterations, [n = 0, dt] () mutable { return dt * n++; });

  array < float > R = zero_order < float, iterations >(x, y0, alpha);

  return 0;
}
