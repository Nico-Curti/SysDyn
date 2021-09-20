// g++ fist_order_kinetic.cpp -std=c++14 -O3 -o fist_order_kinetic

#include <memory>
#include <algorithm>
#include <type_traits>
#include <cassert>

template < class type >
using is_floating_point = typename std :: enable_if < std :: is_floating_point < type > :: value > :: type *;

template < class type, is_floating_point < type > = nullptr >
using array = std :: unique_ptr < type[] >;

/**
* @brief 1st order kinetic
*
* @param x List of time points.
* @param p0 Initial condition of the product.
* @param r0 Initial condition of the reagent.
* @param kf Constant of the forward reaction.
* @param kb Constant of the backward reaction.
* @param P The resulting product array.
* @param R The resulting reagent array.
*
* @tparam type Data-type of arrays
* @tparam Length of time points.
*
*/
template < class type, int32_t N >
void first_order (const array < type > & x, const type & p0, const type & r0,
                  const type & kf, const type & kb,
                  array < type > & P, array < type > & R
                 )
{
  // determine the interval as diff
  const type dx = x[1] - x[0]; // Note: we are assuming it is constant!!

  // Set the initial condition
  P[0] = p0;
  R[0] = r0;

  // Integrate the equation using the Euler method
  for (int32_t i = 0; i < N; ++i)
  {
    P[i + 1] = P[i] + ( kf * R[i] - kb * P[i]) * dx;
    R[i + 1] = R[i] + (-kf * R[i] + kb * P[i]) * dx;
  }
}

int32_t main (/*int32_t argc, char ** argv*/)
{
  const float p0 = 0.f;
  const float r0 = 1.f;
  const float kf = 0.3f;
  const float kb = 0.6f;

  const int32_t iterations = 1000;
  const float dt = 1e-2f;

  array < float > time(new float[iterations]);
  std :: generate_n(time.get(), iterations, [n = 0, dt] () mutable { return dt * n++; });

  array < float > resulting_product(new float[iterations + 1]); // +1 -> initial condition
  array < float > resulting_reagent(new float[iterations + 1]); // +1 -> initial condition
  array < float > resulting_total(new float[iterations + 1]); // +1 -> initial condition

  first_order < float, iterations >(time, p0, r0, kf, kb, resulting_product, resulting_reagent);

  std :: transform (resulting_product.get(), resulting_product.get() + iterations + 1,
                    resulting_reagent.get(),
                    resulting_total.get(), std :: plus < float >());

  assert (std :: all_of(resulting_total.get(), resulting_total.get() + iterations,
                        [] (const float & tot) { return std :: abs(tot - 1.) < 1e-5f; }));

  return 0;
}