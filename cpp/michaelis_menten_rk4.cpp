// g++ michaelis_menten_rk4.cpp -std=c++14 -O3 -o michaelis_menten_rk4

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
* @brief Michaelis Menten kinetic
*
* @param x List of time points.
* @param s0 Initial condition of the substrate.
* @param e0 Initial condition of the enzyme.
* @param es0 Initial condition of the substrate+enzyme.
* @param p0 Initial condition of the product.
* @param kf Constant of the forward reaction.
* @param kb Constant of the backward reaction.
* @param kc Constant of the product reaction.
* @param S The resulting substrate array.
* @param E The resulting enzyme array.
* @param ES The resulting enzyme+substrate array.
* @param P The resulting product array.
*
* @tparam type Data-type of arrays
* @tparam Length of time points.
*
*/
template < class type, int32_t N >
void MichaelisMenten (const array < type > & x,
                      const type & s0, const type & e0, const type & es0, const type & p0,
                      const type & kf, const type & kb, const type &kc,
                      array < type > & S, array < type > & E, array < type > & ES, array < type > & P
                     )
{
  // determine the interval as diff
  const type dx = x[1] - x[0]; // Note: we are assuming it is constant!!

  // Set the initial condition
  S[0] = s0;
  E[0] = e0;
  ES[0] = es0;
  P[0] = p0;

  // set the equation functions
  auto dS = [](const type & S, const type & E, const type & ES,
               const type & kf, const type & kb)
            {
              return -kf * S * E + kb * ES;
            };
  auto dE = [](const type & S, const type & E, const type & ES,
               const type & kf, const type & kb, const type & kc)
            {
              return -kf * S * E + kb * ES + kc * ES;
            };
  auto dES = [](const type & S, const type & E, const type & ES,
                const type & kf, const type & kb, const type & kc)
             {
               return kf * S * E - kb * ES - kc * ES;
             };
  auto dP = [](const type & ES, const type & kc)
            {
              return kc * ES;
            };


  // Integrate the equation using the RK4 method
  for (int32_t i = 0; i < N - 1; ++i)
  {
    const type ks1  = dx * dS( S[i], E[i], ES[i], kf, kb);
    const type ke1  = dx * dE( S[i], E[i], ES[i], kf, kb, kc);
    const type kes1 = dx * dES(S[i], E[i], ES[i], kf, kb, kc);
    const type kp1  = dx * dP(ES[i], kc);

    const type ks2  = dx * dS( S[i] + .5 * ks1, E[i] + .5 * ke1, ES[i] + .5 * kes1, kf, kb);
    const type ke2  = dx * dE( S[i] + .5 * ks1, E[i] + .5 * ke1, ES[i] + .5 * kes1, kf, kb, kc);
    const type kes2 = dx * dES(S[i] + .5 * ks1, E[i] + .5 * ke1, ES[i] + .5 * kes1, kf, kb, kc);
    const type kp2  = dx * dP(ES[i] + .5 * kes1, kc);

    const type ks3  = dx * dS( S[i] + .5 * ks2, E[i] + .5 * ke2, ES[i] + .5 * kes2, kf, kb);
    const type ke3  = dx * dE( S[i] + .5 * ks2, E[i] + .5 * ke2, ES[i] + .5 * kes2, kf, kb, kc);
    const type kes3 = dx * dES(S[i] + .5 * ks2, E[i] + .5 * ke2, ES[i] + .5 * kes2, kf, kb, kc);
    const type kp3  = dx * dP(ES[i] + .5 * kes2, kc);

    const type ks4  = dx * dS( S[i] + ks3, E[i] + ke3, ES[i] + kes3, kf, kb);
    const type ke4  = dx * dE( S[i] + ks3, E[i] + ke3, ES[i] + kes3, kf, kb, kc);
    const type kes4 = dx * dES(S[i] + ks3, E[i] + ke3, ES[i] + kes3, kf, kb, kc);
    const type kp4  = dx * dP(ES[i] + kes3, kc);

    S[i + 1]  = S[i]  + type(1. / 6.) * (ks1  + type(2.) * ks2  + type(2.) * ks3  + ks4);
    E[i + 1]  = E[i]  + type(1. / 6.) * (ke1  + type(2.) * ke2  + type(2.) * ke3  + ke4);
    ES[i + 1] = ES[i] + type(1. / 6.) * (kes1 + type(2.) * kes2 + type(2.) * kes3 + kes4);
    P[i + 1]  = P[i]  + type(1. / 6.) * (kp1  + type(2.) * kp2  + type(2.) * kp3  + kp4);
  }
}


int32_t main (/*int32_t argc, char ** argv*/)
{
  const float s0  = 10.f;
  const float e0  = 1.f;
  const float es0 = 0.f;
  const float p0  = 0.f;

  const float kf = 1.f;
  const float kb = 1e-2f;
  const float kc = 1.f;

  const int32_t iterations = 101;
  const float dt = (20.f - 1e-2f) / iterations;

  array < float > time(new float[iterations]);
  std :: generate_n(time.get(), iterations, [n = 0, dt] () mutable { return 1e-2f + dt * n++; });

  array < float > S(new float[iterations]);
  array < float > E(new float[iterations]);
  array < float > ES(new float[iterations]);
  array < float > P(new float[iterations]);

  MichaelisMenten < float, iterations >(time, s0, e0, es0, p0,
                                        kf, kb, kc,
                                        S, E, ES, P);

  return 0;
}