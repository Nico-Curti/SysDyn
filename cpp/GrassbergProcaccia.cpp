#include <iostream>
#include <algorithm>
#include <climits>
#include <cmath>
#define LOG2 0.69314718055994529

template<typename T> T* GrassbergerProcaccia(T *x, T *y, const int &N, int omit_pts = 3)
{
	T min_eps =  std::numeric_limits<T>::infinity(),
	  max_eps = -std::numeric_limits<T>::infinity(),
	  *ED = new T[N*(N-1)/2],
	  *eps_vec = nullptr,
	  Npairs,
	  *C_eps = nullptr,
	  sx = (T)0., 
	  sy = (T)0., 
	  sxy = (T)0., 
	  sx2 = (T)0., 
	  w,
	  xp, yp,
	  *parameters = new T[4];
	int eps_vec_size,
		k1 = omit_pts,
		k2,
		idx;

#pragma omp parallel for reduction (min: min_eps) reduction(max: max_eps) private(idx)
	for(int i = 0; i < N; ++i)
		for(int j = i+1; j < N; ++j)
		{
			idx = (N * (N - 1) / 2) - (N - i) * ((N - i) - 1) / 2 + j - i - 1;
			ED[idx] = (T)std::sqrt( (x[i] - x[j])*(x[i] - x[j]) + (y[i] - y[j])*(y[i] - y[j]) );
			min_eps = (ED[idx] < min_eps && ED[idx] != 0.f) ? ED[idx] : min_eps;
			max_eps = (ED[idx] > max_eps && ED[idx] != 0.f) ? ED[idx] : max_eps;
		}
	max_eps = (T)std::pow(2, std::ceil(std::log(max_eps) / LOG2));	
	eps_vec_size = (int)std::floor( (std::log(max_eps/min_eps) / LOG2) ) + 1;
	eps_vec = new T[eps_vec_size];
	std::generate(eps_vec, eps_vec + eps_vec_size, [n = 0, &max_eps]() mutable {return max_eps*(T)std::pow(2, - n++);});
	Npairs = N * (N - 1) / 2;

	C_eps = new T[eps_vec_size];
	std::transform(eps_vec, eps_vec + eps_vec_size,
				   C_eps,
				   [&Npairs, &ED, &N](const T &eps)
				   {
				   	return (T)2.*std::count_if(ED, ED + Npairs,
				   							[&eps](const T &ed)
				   							{
				   								return (ed < eps) ? 1. : 0.;
				   							}) / Npairs;
				   });
	k2 = eps_vec_size - omit_pts;
	// Compute correlation dimension as linear fit (slope)
#pragma omp parallel for reduction(+ : sx, sy, sxy, sx2)
    for(int i = k1; i < k2; ++i)
    {
    	xp = (T)std::log(eps_vec[i]) / LOG2;
    	yp = (T)std::log(C_eps[i]) / LOG2;
        sx += xp;
        sy += yp;
        sxy += xp * yp;
        sx2 += xp * xp;
    }
    w = (k2 - k1) * sx2 - sx * sx;
    parameters[0] = ((k2 - k1) * sxy - sx*sy) / w; //slope
    parameters[1] = (T)std::sqrt((k2 - k1) / w); // err slope
    parameters[2] = (sx2 * sy - sxy*sx) / w; //intercept
    parameters[3] = (T)std::sqrt(sx2 / w); // err intercept

	delete[] ED;
	delete[] eps_vec;
	delete[] C_eps;

    return parameters;
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

	double *params = GrassbergerProcaccia(x, y, Npts);
	std::cout << "slope : " << params[0] << std::endl
			  << "err slope : " << params[1] << std::endl
			  << "intercept : " << params[2] << std::endl
			  << "err intercept : " << params[3] << std::endl
			  << std::endl;

	delete[] x;
	delete[] y;

	return 0;
}
