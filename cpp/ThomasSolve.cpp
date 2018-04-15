#include <iostream>
#include <iterator>

template<typename T> T* Thomas(const T* b, const T* a, T *c, T *d, const int &nb)
{
    int i = 0;
    T *x = new T[nb], id;
    c[i] /= b[i];
    d[i] /= b[i];
    for(i = 1; i < nb-1; ++i)
    {
        id = (T)1. / (b[i] - a[i-1] * c[i-1]);
        d[i] = (d[i] - a[i-1] * d[i-1]) * id;
        c[i] *= id;
    }
    d[i] = (d[i] - a[i-1] * d[i-1]) / (b[i] - a[i-1] * c[i-1]);
    x[nb - 1] = d[i];
    for(i = nb-2; i != -1; --i) x[i] = d[i] - c[i] * x[i+1];
    return x;
}

template<typename T> T* Thomas(const T* b, T* c, T* d, const int &nb)
{
    int i = 0;
    T *x = new T[nb], id, a = c[i];
    c[i] /= b[i];
    d[i] /= b[i];
    for(i = 1; i < nb-1; ++i)
    {
        id = (T)1. / (b[i] - a*c[i-1]);
        d[i] = (d[i] - a*d[i-1]) * id;
        a = c[i];
        c[i] *= id;
    }
    d[i] = (d[i] - a*d[i-1]) / (b[i] - a * c[i-1]);
    x[nb - 1] = d[i];
    for(i = nb-2; i != -1; --i) x[i] = d[i] - c[i] * x[i+1];
    return x;
}


int main(int argc, char **argv)
{
    float *a = new float[2],
          *b = new float[3],
          *c = new float[2],
          *d = new float[3];
    a[0] = 4.f; a[1] = 3.f;
    b[0] = 9.f; b[1] = -7.f; b[2]  = 8.f;
    c[0] = 1.f; c[1] = 2.f;
    d[0] = 5.f; d[1] = 6.f; d[2] = 2.f;

    float *x = Thomas(b, a, c, d, 3);
    std::cout << "Solution:" << std::endl;
    std::copy(x, x + 3, std::ostream_iterator<float>(std::cout, " "));
    std::cout << std::endl;

    delete[] a;
    delete[] b;
    delete[] c;
    delete[] d;
    delete[] x;
    return 0;
}