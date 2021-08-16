#include <iostream>
#include <iterator>
#include <memory>
#include <algorithm>

std :: unique_ptr < float[] > Thomas (const float * b, const float * a,
                                      float * c, float * d,
                                      const int32_t & nb)
{
  int32_t i = 0;
  std :: unique_ptr < float[] > x(new float[nb]);

  float tmp = c[i];

  c[i] /= b[i];
  d[i] /= b[i];

  for (i = 1; i < nb - 1; ++i)
  {
    const float scale = a ? a[i - 1] : tmp;
    const float id = 1.f / (b[i] - scale * c[i-1]);
    d[i] = (d[i] - a[i-1] * d[i-1]) * id;
    tmp = c[i];
    c[i] *= id;
  }

  tmp = a ? a[i - 1] : tmp;
  d[i] = (d[i] - tmp * d[i-1]) / (b[i] - tmp * c[i-1]);
  x[nb - 1] = d[i];

  for (i = nb-2; i != -1; --i)
    x[i] = d[i] - c[i] * x[i+1];

  return x;
}


int main (int argc, char ** argv)
{
  std :: array < float, 2 > a = {4.f, 3.f};
  std :: array < float, 3 > b = {9.f, -7.f, 8.f};
  std :: array < float, 2 > c = {1.f, 2.f};
  std :: array < float, 3 > d = {5.f, 6.f, 2.f};

  auto x = Thomas(b.data(), a.data(), c.data(), d.data(), 3);

  std :: cout << "Solution:" << std :: endl;
  std :: copy_n(x.get(), 3, std :: ostream_iterator < float >(std :: cout, " "));
  std :: cout << std :: endl;

  return 0;
}
