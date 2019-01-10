#include <iostream>
#include <numeric>
#include <algorithm>
#include <iterator>
#include <array>
#include <memory>
#include <cmath>

auto GaussSeidel(const std::array<std::array<float, 4>, 4> &A, const std::array<float, 4> &b, const int &N, const int &M)
{
  const int ITERATION_LIMIT = static_cast<int>(1e20);
  int cnt = 0;
  std::unique_ptr<float[]> x_new(new float[M]),
                           res(new float[M]);
  float diff;
  for(int it = 0; it < ITERATION_LIMIT; ++it)
  {
#ifdef DEBUG
    std::cout << "Current solution: ";
    std::copy_n(x_new.get(), M, std::ostream_iterator<float>(std::cout, " "));
    std::cout << std::endl;
#endif
    cnt = 0;
    for(int i = 0; i < N; ++i)
    {
      diff = std::inner_product(A[i].begin(), A[i].begin() + i, x_new.get(), 0.f) +
             std::inner_product(A[i].begin() + i + 1, A[i].begin() + M, x_new.get() + i + 1, 0.f);
      x_new[i] = (b[i] - diff) / A[i][i];
      cnt += (std::fabs(res[i] - x_new[i]) > 1e-8) ? 1 : 0;
    }
    if(!cnt) break;
    std::copy_n(res.get(), M, x_new.get());
  }
  return res;
}


int main(int argc, char **argv)
{
  std::array<std::array<float, 4>, 4> A;
  std::array<float, 4> b,
                       errors;

  A[0][0] = 10.f; A[0][1] = -1.f; A[0][2] = 2.f;  A[0][3] = 0.f;
  A[1][0] = -1.f; A[1][1] = 11.f; A[1][2] = -1.f; A[1][3] = 3.f;
  A[2][0] = 2.f;  A[2][1] = -1.f; A[2][2] = 10.f; A[2][3] = -1.f;
  A[3][0] = 0.f;  A[3][1] = 3.f;  A[3][2] = -1.f; A[3][3] = 8.f;

  b[0] = 6.f; b[1] = 25.f; b[2] = -11.f; b[3] = 15.f;

  std::cout << "System:" << std::endl;
  for(int i = 0; i < 4; ++i)
  {
    for(int j = 0; j < 3; ++j) std::cout << A[i][j] << "*x" << j << " + ";
    std::cout << A[i][3] << "*x" << 3 << " = " << b[i] << std::endl;
  }

  auto res = GaussSeidel(A, b, 4, 4);
  std::cout << "Solution:" << std::endl;
  std::copy_n(res.get(), 4, std::ostream_iterator<float>(std::cout, " "));
  std::cout << std::endl;
  std::transform(A.begin(), A.end(),
                 b.begin(), errors.begin(),
                 [&](const std::array<float, 4> &Ai, const float &bi)
                 {
                   return std::inner_product(Ai.begin(), Ai.end(), res.get(), 0.f) - bi;
                 });
  std::cout << "Errors:" << std::endl;
  std::copy_n(errors.begin(), 4, std::ostream_iterator<float>(std::cout, " "));
  std::cout << std::endl;

  return 0;
}
