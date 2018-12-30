#include <iostream> // std::cout
#include <cstring> // std::memset
#include <random> // std::uniform_real_distribution
#include <utility> // std::pair
#include <algorithm>
#include <numeric>
#include <array>
#ifdef STDOUT
#include <fstream>
#endif

static constexpr int N = 20;
static constexpr int SEED = 123;

std::mt19937 eng(SEED);
std::uniform_real_distribution<float> distr(0.f, 1.f);

auto kalman_filter = [&](float *x, float *P,
                         const float *F, const float *H,
                         const float *motion,
                         const float *Q, const float &R,
                         const float &obs_x, const float &obs_y) -> std::pair<float, float>
  {
    // Distance between measured and current position-belief
    // y = measurements.T - H * x
    float y0 = obs_x - (H[0]*x[0] + H[1]*x[1] + H[2]*x[2] + H[3]*x[3]),
          y1 = obs_y - (H[4]*x[0] + H[5]*x[1] + H[6]*x[2] + H[7]*x[3]);
    // P * H.T
    float phT00 = P[0]* H[0] + P[1]* H[1] + P[2] *H[2] + P[3] *H[3],
          phT01 = P[0]* H[4] + P[1]* H[5] + P[2] *H[6] + P[3] *H[7],
          phT10 = P[4]* H[0] + P[5]* H[1] + P[6] *H[2] + P[7] *H[3],
          phT11 = P[4]* H[4] + P[5]* H[5] + P[6] *H[6] + P[7] *H[7],
          phT20 = P[8]* H[0] + P[9]* H[1] + P[10]*H[2] + P[11]*H[3],
          phT21 = P[8]* H[4] + P[9]* H[5] + P[10]*H[6] + P[11]*H[7],
          phT30 = P[12]*H[0] + P[13]*H[1] + P[14]*H[2] + P[15]*H[3],
          phT31 = P[12]*H[4] + P[13]*H[5] + P[14]*H[6] + P[15]*H[7];
    // S = H * P * H.T + R  // Residual covariance
    float S00  = H[0] * phT00 + H[1] * phT10 + H[2] * phT20 + H[3] * phT30 + R,
          S01  = H[0] * phT01 + H[1] * phT11 + H[2] * phT21 + H[3] * phT31 + R,
          S10  = H[4] * phT00 + H[5] * phT10 + H[6] * phT20 + H[7] * phT30 + R,
          S11  = H[4] * phT01 + H[5] * phT11 + H[6] * phT21 + H[7] * phT31 + R;
    float inverse_scale = 1.f / (S00*S11 - S01*S10);
    // K = P * H.T * S.inv  // Kalman gain
    float K00 = inverse_scale * (phT00 *  S11 + phT01 * -S10),
          K01 = inverse_scale * (phT00 * -S01 + phT01 *  S00),
          K10 = inverse_scale * (phT10 *  S11 + phT11 * -S10),
          K11 = inverse_scale * (phT10 * -S01 + phT11 *  S00),
          K20 = inverse_scale * (phT20 *  S11 + phT21 * -S10),
          K21 = inverse_scale * (phT20 * -S01 + phT21 *  S00),
          K30 = inverse_scale * (phT30 *  S11 + phT31 * -S10),
          K31 = inverse_scale * (phT30 * -S01 + phT31 *  S00);
    // I - K*H
    float KH00 = 1.f - (K00 * H[0] + K01 * H[4]),
          KH01 =     - (K00 * H[1] + K01 * H[5]),
          KH02 =     - (K00 * H[2] + K01 * H[6]),
          KH03 =     - (K00 * H[3] + K01 * H[7]),
          KH10 =     - (K10 * H[0] + K11 * H[4]),
          KH11 = 1.f - (K10 * H[1] + K11 * H[5]),
          KH12 =     - (K10 * H[2] + K11 * H[6]),
          KH13 =     - (K10 * H[3] + K11 * H[7]),
          KH20 =     - (K20 * H[0] + K21 * H[4]),
          KH21 =     - (K20 * H[1] + K21 * H[5]),
          KH22 = 1.f - (K20 * H[2] + K21 * H[6]),
          KH23 =     - (K20 * H[3] + K21 * H[7]),
          KH30 =     - (K30 * H[0] + K31 * H[4]),
          KH31 =     - (K30 * H[1] + K31 * H[5]),
          KH32 =     - (K30 * H[2] + K31 * H[6]),
          KH33 = 1.f - (K30 * H[3] + K31 * H[7]);
  // P = (I - K*H)*P
  float Pnew00 = KH00*P[0] + KH01*P[4] + KH02*P[8]  + KH03*P[12],
        Pnew01 = KH00*P[1] + KH01*P[5] + KH02*P[9]  + KH03*P[13],
        Pnew02 = KH00*P[2] + KH01*P[6] + KH02*P[10] + KH03*P[14],
        Pnew03 = KH00*P[3] + KH01*P[7] + KH02*P[11] + KH03*P[15],
        Pnew10 = KH10*P[0] + KH11*P[4] + KH12*P[8]  + KH13*P[12],
        Pnew11 = KH10*P[1] + KH11*P[5] + KH12*P[9]  + KH13*P[13],
        Pnew12 = KH10*P[2] + KH11*P[6] + KH12*P[10] + KH13*P[14],
        Pnew13 = KH10*P[3] + KH11*P[7] + KH12*P[11] + KH13*P[15],
        Pnew20 = KH20*P[0] + KH21*P[4] + KH22*P[8]  + KH23*P[12],
        Pnew21 = KH20*P[1] + KH21*P[5] + KH22*P[9]  + KH23*P[13],
        Pnew22 = KH20*P[2] + KH21*P[6] + KH22*P[10] + KH23*P[14],
        Pnew23 = KH20*P[3] + KH21*P[7] + KH22*P[11] + KH23*P[15],
        Pnew30 = KH30*P[0] + KH31*P[4] + KH32*P[8]  + KH33*P[12],
        Pnew31 = KH30*P[1] + KH31*P[5] + KH32*P[9]  + KH33*P[13],
        Pnew32 = KH30*P[2] + KH31*P[6] + KH32*P[10] + KH33*P[14],
        Pnew33 = KH30*P[3] + KH31*P[7] + KH32*P[11] + KH33*P[15];
  float FP00 = F[0]*Pnew00 + F[1]*Pnew10 + F[2]*Pnew20 + F[3]*Pnew30,
        FP01 = F[0]*Pnew01 + F[1]*Pnew11 + F[2]*Pnew21 + F[3]*Pnew31,
        FP02 = F[0]*Pnew02 + F[1]*Pnew12 + F[2]*Pnew22 + F[3]*Pnew32,
        FP03 = F[0]*Pnew03 + F[1]*Pnew13 + F[2]*Pnew23 + F[3]*Pnew33,
        FP10 = F[4]*Pnew00 + F[5]*Pnew10 + F[6]*Pnew20 + F[7]*Pnew30,
        FP11 = F[4]*Pnew01 + F[5]*Pnew11 + F[6]*Pnew21 + F[7]*Pnew31,
        FP12 = F[4]*Pnew02 + F[5]*Pnew12 + F[6]*Pnew22 + F[7]*Pnew32,
        FP13 = F[4]*Pnew03 + F[5]*Pnew13 + F[6]*Pnew23 + F[7]*Pnew33,
        FP20 = F[8]*Pnew00 + F[9]*Pnew10 + F[10]*Pnew20 + F[11]*Pnew30,
        FP21 = F[8]*Pnew01 + F[9]*Pnew11 + F[10]*Pnew21 + F[11]*Pnew31,
        FP22 = F[8]*Pnew02 + F[9]*Pnew12 + F[10]*Pnew22 + F[11]*Pnew32,
        FP23 = F[8]*Pnew03 + F[9]*Pnew13 + F[10]*Pnew23 + F[11]*Pnew33,
        FP30 = F[12]*Pnew00 + F[13]*Pnew10 + F[14]*Pnew20 + F[15]*Pnew30,
        FP31 = F[12]*Pnew01 + F[13]*Pnew11 + F[14]*Pnew21 + F[15]*Pnew31,
        FP32 = F[12]*Pnew02 + F[13]*Pnew12 + F[14]*Pnew22 + F[15]*Pnew32,
        FP33 = F[12]*Pnew03 + F[13]*Pnew13 + F[14]*Pnew23 + F[15]*Pnew33;
    // Update x
    x[0] += K00*y0 + K01*y1; x[1] += K10*y0 + K11*y1; x[2] += K20*y0 + K21*y1; x[3] += K30*y0 + K31*y1;

    x[0] = F[0] *x[0] + F[1] *x[1] + F[2] *x[2] + F[3] *x[3] + motion[0];
    x[1] = F[4] *x[0] + F[5] *x[1] + F[6] *x[2] + F[7] *x[3] + motion[1];
    x[2] = F[8] *x[0] + F[9] *x[1] + F[10]*x[2] + F[11]*x[3] + motion[2];
    x[3] = F[12]*x[0] + F[13]*x[1] + F[14]*x[2] + F[15]*x[3] + motion[3];

    P[0] = FP00*F[0]  + FP01*F[1]  + FP02*F[2]  + FP03*F[3]  + Q[0];
    P[1] = FP00*F[4]  + FP01*F[5]  + FP02*F[6]  + FP03*F[7]  + Q[1];
    P[2] = FP00*F[8]  + FP01*F[9]  + FP02*F[10] + FP03*F[11] + Q[2];
    P[3] = FP00*F[12] + FP01*F[13] + FP02*F[14] + FP03*F[15] + Q[3];

    P[4] = FP10*F[0] + FP11 *F[1]  + FP12*F[2]  + FP13*F[3]  + Q[4];
    P[5] = FP10*F[4] + FP11 *F[5]  + FP12*F[6]  + FP13*F[7]  + Q[5];
    P[6] = FP10*F[8] + FP11 *F[9]  + FP12*F[10] + FP13*F[11] + Q[6];
    P[7] = FP10*F[12] + FP11*F[13] + FP12*F[14] + FP13*F[15] + Q[7];

    P[8]  = FP20*F[0] + FP21 *F[1]  + FP22*F[2]  + FP23*F[3]  + Q[8];
    P[9]  = FP20*F[4] + FP21 *F[5]  + FP22*F[6]  + FP23*F[7]  + Q[9];
    P[10] = FP20*F[8] + FP21 *F[9]  + FP22*F[10] + FP23*F[11] + Q[10];
    P[11] = FP20*F[12] + FP21*F[13] + FP22*F[14] + FP23*F[15] + Q[11];

    P[12] = FP30*F[0] + FP31 *F[1]  + FP32*F[2]  + FP33*F[3]  + Q[12];
    P[13] = FP30*F[4] + FP31 *F[5]  + FP32*F[6]  + FP33*F[7]  + Q[13];
    P[14] = FP30*F[8] + FP31 *F[9]  + FP32*F[10] + FP33*F[11] + Q[14];
    P[15] = FP30*F[12] + FP31*F[13] + FP32*F[14] + FP33*F[15] + Q[15];

    return std::make_pair(x[0], x[1]);
  };


int main(int argc, char **argv)
{
  std::array<float, 4> initial, // initial state 4-tuple of location and velocity: (x0, x1, x0_dot, x1_dot)
                       motion;  // external motion added to state vector x
  std::array<float, 16> P,      // initial uncertainty
                        Q;      // motion noise (same shape as P)

  std::array<float, N> true_x,
                       true_y,
                       obs_x,
                       obs_y;
  std::array<float, 8> H,       // measurement function: position = H*x
                       F;       // next state function: x_prime = F*x

  constexpr float R = 1e-1f * 1e-1f; // measurement noise

  std::fill_n(initial.begin(), 4,  0.f);
  std::fill_n(motion.begin(), 4,  0.f);

  std::fill_n(P.begin(), 16, 0.f); // identity *1000
  P[0] = 1e3f; P[5] = 1e3f; P[10] = 1e3f; P[15] = 1e3f;

  std::memset(Q.begin(), 16, 0.f); // Identity
  Q[0] = 1.f; Q[5] = 1.f; Q[10] = 1.f; Q[15] = 1.f;

  for(int i = 0; i < N; ++i)
  {
    true_x[i] = static_cast<float>(i) / N;
    true_y[i] = true_x[i]*true_x[i];
    obs_x[i]  = true_x[i] + .05f*distr(eng)*true_x[i];
    obs_y[i]  = true_y[i] + .05f*distr(eng)*true_y[i];
    //obs_x[i]  = true_x[i] + .05f*((float)std::rand() / RAND_MAX)*true_x[i];
    //obs_y[i]  = true_y[i] + .05f*((float)std::rand() / RAND_MAX)*true_y[i];
  }

  F[0]  = 1.f; F[1]  = 0.f; F[2]  = 1.f; F[3]  = 0.f;
  F[4]  = 0.f; F[5]  = 1.f; F[6]  = 0.f; F[7]  = 1.f;
  F[8]  = 0.f; F[9]  = 0.f; F[10] = 1.f; F[11] = 0.f;
  F[12] = 0.f; F[13] = 0.f; F[14] = 0.f; F[15] = 1.f;

  H[0] = 1.f; H[1] = 0.f; H[2] = 0.f; H[3] = 0.f;
  H[4] = 0.f; H[5] = 1.f; H[6] = 0.f; H[7] = 0.f;

  std::array<std::pair<float, float>, N> rec;
  std::transform(obs_x.begin(), obs_x.end(),
                 obs_y.begin(), rec.begin(),
                 [&](const float &obx, const float &oby)
                 {
                  return kalman_filter(initial.data(), P.data(), F.data(), H.data(), motion.data(), Q.data(), R.data(), obx.data(), oby.data());
                 });

  std::cout << "\tTrue coord\tFiltered coord" << std::endl;
  for(int i = 0; i < N; ++i)
    std::cout << "(x, y) " << obs_x[i] << ", " << obs_y[i] << " -> " << rec[i].first << ", " << rec[i].second << std::endl;

#ifdef STDOUT
  std::ofstream os("kalman_points.dat");
  os << "obs_x,obs_y,rec_x,rec_y" << std::endl;
  for(int i = 0; i < N; ++i)
    os << obs_x[i] << "," << obs_y[i] << "," << rec[i].first << "," << rec[i].second << std::endl;
  os.close();

  os.open("kalman_points.py")
  os << "#!/usr/bin/python" << std::endl
     << "import pandas as pd" << std::endl
     << "import matplotlib.pylab as plt" << std::endl << std::endl
     << "pts = pd.read_csv('kalman_points.dat', sep=',', header=0)" << std::endl
     << "fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8,8))" << std::endl
     << "ax.plot(pts.obs_x, pts.obs_y, 'bo', lw=1, label='observed signal')" << std::endl
     << "ax.plot(pts.rec_x, pts.rec_y, 'r-', lw=1, label='kalman signal')" << std::endl
     << "plt.legend(loc='best', fontsize=14)" << std::endl
     << "plt.savefig('kalman_points.png')" << std::endl;
  os.close();
  std::system("python ./kalman_points.py")
#endif
  return 0;
}
