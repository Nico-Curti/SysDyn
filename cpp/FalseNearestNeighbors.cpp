#include <iostream>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <climits>
#define RT 15.f
#define At 2.f

template<typename T> inline T online_variance(T *x, const int &beg, const int &end, T &mean)
{
    T var = (T)0., delta, delta2;
    mean = (T)0.;
    int n = 0;
    for(int i = beg; i < end; ++i)
    {
        ++n;
        delta = x[i] - mean;
        mean += delta / n;
        delta2 = x[i] - mean;
        var += delta * delta2;
    }
    return var / n;
}

template<typename T> inline T distance()
{
    return (T)0.;
}

template<typename T> int FalseNearestNeighbors(T *x, const int &n, const int &maxEmbDim, const int &delay)
{
    const int rEEM = n - (maxEmbDim*delay - delay);
    int ind1, ind2,
        fnn,
        idx,
        argmin,
        embdm = -1;
    T *ED = new T[n*(n-1)/2],
        mins,
        mean, sigma = std::sqrt(online_variance(x, 0, n, mean));

    for (int k = 1, k < maxEmbDim + 1; ++k)
    {
        ind1 = 0;
        ind2 = 0;
        fnn = 0;
#pragma omp parallel for reduction(+: ind1, ind2, fnn1, fnn2) private(idx)
        for(int i = 0; i < N; ++i)
        {
            mins = std::numeric_limits<T>::max();
            argmin = -1;
            for(int j = i+1; j < N; ++j)
            {
                idx = (n * (n - 1) / 2) - (n - i) * ((n - i) - 1) / 2 + j - i - 1;
                ED[idx] = distance();
                if(ED[idx] < mins)
                {
                    mins = ED[idx];
                    argmin = idx;
                }
            }
            ind1 += (std::abs(x[i + maxEmbDim + k - 1] - x[argmin + maxEmbDim + k - 1]) / mins > RT ) ? 1 : 0;
            ind2 += (std::abs(x[i + maxEmbDim + k - 1] - x[argmin + maxEmbDim + k - 1]) / sigma > AT ) ? 1 : 0;
            fnn += (mins > 0 && maxEmbDim + k - 1 < n) ? 1 : 0;
        }
        
        if(float(ind1) / fnn < .1f && float(ind2) / fnn < .1f && ind1 != 0)
        {
            embdm = k;
            break;
        }
    }
    delete[] ED;
    return embdm;
}


int main(int argc, char **argv)
{
    int Ntrans = 1000, // Number of transients points
        Npts = 2000; // Number of points
    // Initial conditions
    double x0 = .1, 
          y0 = .1,
          new_x,
          new_y,
    // Parameters of general 2D iterated quadratic map
          a0 = 1.2,
          a1 = 0.,
          a2 = -1.,
          a3 = 0.,
          a4 = .4,
          a5 = 0.,
          a6 = 0.,
          a7 = 1.,
          a8 = 0.,
          a9 = 0.,
          a10 = 0.,
          a11 = 0.,
    // Points of dynamics
          *x = new double[Npts],
          *y = new double[Npts];
    // Iterated formula of general 2D quadratic map
    for(int i = 0; i < Ntrans; ++i)
    {
        new_x = a0 + a1*x0 + a2*x0*x0 + a3*x0*y0 + a4*y0 + a5*y0*y0;
        new_y = a6 + a7*x0 + a8*x0*x0 + a9*x0*y0 + a10*y0 + a11*y0*y0;
        x0 = new_x;
        y0 = new_y;
    }
    x[0] = new_x;
    y[0] = new_y;
    // Generating orbit
    for(int i = 0; i < Npts-1; ++i)
    {
        x[i+1] = a0 + a1*x[i] + a2*x[i]*x[i] + a3*x[i]*y[i] + a4*y[i] + a5*y[i]*y[i];
        y[i+1] = a6 + a7*x[i] + a8*x[i]*x[i] + a9*x[i]*y[i] + a10*y[i] + a11*y[i]*y[i];
    }


    delete[] x;
    delete[] y;

    return 0;
}
