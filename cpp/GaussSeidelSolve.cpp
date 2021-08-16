#include <iostream>
#include <numeric>
#include <algorithm>
#include <iterator>
#include <array>
#include <memory>
#include <cmath>

std :: unique_ptr < float[] > GaussSeidel (const std :: array < std :: array < float, 4 >, 4 > & A,
                                           const std :: array < float, 4 > & b,
                                           const int32_t & N, const int32_t & M)
{
  const int32_t ITERATION_LIMIT = static_cast < int32_t >(1e20);
  int32_t cnt = 0;

  std :: unique_ptr < float[] > x_new(new float[M]);
  std :: unique_ptr < float[] > res(new float[M]);

  std :: fill_n(res.get(), M, 0.f);

  float diff;

  for (int32_t it = 0; it < ITERATION_LIMIT; ++it)
  {
    cnt = 0;

    for (int32_t i = 0; i < N; ++i)
    {
      diff = std :: inner_product(A[i].begin(), A[i].begin() + i, x_new.get(), 0.f) +
             std :: inner_product(A[i].begin() + i + 1, A[i].begin() + M, x_new.get() + i + 1, 0.f);
      x_new[i] = (b[i] - diff) / A[i][i];
      cnt += (std::fabs(res[i] - x_new[i]) > 1e-8) ? 1 : 0;
    }

    if ( ! cnt )
      break;

    std :: copy_n(x_new.get(), M, res.get());

    std :: cout << "Current solution: ";
    std :: copy_n(res.get(), M, std :: ostream_iterator < float >(std :: cout, " "));
    std :: cout << std :: endl;
  }

  return res;
}


int main (int argc, char ** argv)
{
  std :: array < std :: array < float, 4 >, 4 > A = {{
                                                      {{10.f, -1.f, 2.f,  0.f}},
                                                      {{-1.f, 11.f, -1.f, 3.f}},
                                                      {{2.f,  -1.f, 10.f, -1.f}},
                                                      {{0.f,  3.f,  -1.f, 8.f}}
                                                    }};

  std :: array < float, 4 > b = {6.f, 25.f, -11.f, 15.f};
  std :: array < float, 4 > errors;

  std :: cout << "System:" << std :: endl;

  for (int32_t i = 0; i < 4; ++i)
  {
    for (int32_t j = 0; j < 3; ++j)
      std :: cout << A[i][j] << "*x" << j << " + ";

    std :: cout << A[i][3] << "*x" << 3 << " = " << b[i] << std :: endl;
  }

  auto res = GaussSeidel(A, b, 4, 4);

  std :: cout << "Solution:" << std :: endl;
  std :: copy_n(res.get(), 4, std :: ostream_iterator < float >(std :: cout, " "));
  std :: cout << std :: endl;
  std :: transform(A.begin(), A.end(),
                   b.begin(), errors.begin(),
                   [&](const std :: array < float, 4 > & Ai, const float & bi)
                   {
                     return std :: inner_product(Ai.begin(), Ai.end(), res.get(), 0.f) - bi;
                   });

  std :: cout << "Errors:" << std :: endl;
  std :: copy_n(errors.begin(), 4, std :: ostream_iterator < float >(std :: cout, " "));
  std :: cout << std :: endl;

  return 0;
}
