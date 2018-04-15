#include <iostream>
#include <numeric>
#include <algorithm>
#include <vector>

#define A 2.f
#define B 5.2f
#define dt .01f
#define ITER 10
#define NUM_REAC 4
/*
                       Brusselator - Stochastic Approach

                                        A -> X
                                    B + X -> Y + C
                                   2X + Y -> 3X
                                        X -> D

We assume that A, B, C and D concetrations are constants, so we can model the
system evolution as a cinetic model reactions

1 reaction :  dX/dt = +V1    dY/dt = 0
2 reaction :  dX/dt = -V2    dY/dt = +V2
3 reaction :  dX/dt = +V3    dY/dt = -V3
4 reaction :  dX/dt = -V4    dY/dt = 0

Considering positive and negative terms of each equation and using the Law of Mass Action
we rewrite the equations as
                        
                                dX/dt = A - BX + X^2Y - X
                                dY/dt = BX - X^2Y

With a stochastic approach, considering a number of molecula about the Avogadro number we can
write a Chemical Master Equation of the Model. Assuming as starting point the state n we
have 4 possible state around that as shown in the scheme:  

            o (x-1 , y+1) 
            
            o (x-1, y)            o (x, y)          o (x+1 , y)
            
                                                    o (x+1, y-1)

Each step will be correlated with a particular coefficient of ricombination (step back) and generation (step forward).
A key parameter of the Gillespie algorithm is the system dimension (omega) which represent the connection
between the deterministic and stochastic model.
In the CME model the omega parameter contributes in the reaction step with 2 or more molecular species.
In this case only the reaction 2 and 3 verify this conditions: in the first case we have to divide the constant
reaction (= 1 in our case) for the omega value (1 specie) while in the second case there are 2 species so 
we divide the value for the square of omega.
In conclusion the CME system will be: 

                                c1 = k1 = 1
                                c2 = k2/Omega = 1/Omega
                                c3 = 2k3/Omega^2 = 2/Omega^2
                                c4 = k4 = 1

                dP(x,y,t)/dt = -(c1A + c2BX + c3X^2Y + c4X)P(x,y,t)
                                + c1A P(x-1,y,t)
                                + c2B(X+1) P(x+1, y-1, t)
                                + c3(X-1)(X-2)(Y+1) P(x-1, y+1, t)
                                + c4(X+1) P(x+1,y, t)

*/

template<typename T> inline void reaction(const T &x, const T &y, T *w, const T &noise)
{
	w[0] = A * noise;
	w[1] = B * x;
	w[2] = x;
	w[3] = x * (x - 1) * y / (noise * noise);
}

template<typename T> inline void reac_step(T &x, T &y, T *comp, const T &act_step)
{
	if(act_step < comp[0]) ++x;
	else if(act_step < comp[1]) {--x; ++y;}
	else if(act_step < comp[2]) --x;
	else if(act_step < comp[3]) {++x; --y;}
}

template<typename T> void CME(std::vector<T> &x, std::vector<T> &y, const T &omega, unsigned int seed = 123)
{
	T trans = (T)0., told = (T)0., t = (T)0., 
	  *w = new T[NUM_REAC],
	  *c = new T[NUM_REAC],
	  z1, z2, tau, uct, tmpx = x[0], tmpy = y[0];
	std::srand(seed);

	while(t <= ITER + trans)
	{
#ifdef DEBUG
		std::cout << "Time " << t << std::endl;
#endif
		reaction(tmpx, tmpy, w, omega);
		c[0] = w[0];
		for(int j = 1; j < NUM_REAC; ++j) c[j] = c[j-1] + w[j];

		z1 = static_cast <T> (std::rand()) / static_cast <T> (RAND_MAX);
		z2 = static_cast <T> (std::rand()) / static_cast <T> (RAND_MAX);
		tau = - std::log(z1) / c[NUM_REAC - 1];
		uct = z2 * c[NUM_REAC - 1];
		t += tau;
		reac_step(tmpx, tmpy, c, uct);
		if(t > trans && t > told + dt)
		{
			x.push_back(tmpx);
			y.push_back(tmpy);
			told = t;
		}
	}
}

#define DEBUG
int main(int argc, char **argv)
{
	float omega = 1000.f;
	std::vector<float> x, y;
	x.push_back(1.6f);
	y.push_back(2.8f);

	CME(x, y, omega);

#ifdef DEBUG
	std::cout << "x = [";
	for(const auto &i : x) std::cout << i << ",";
	std::cout << "]" << std::endl;
	std::cout << "y = [";
	for(const auto &i : y) std::cout << i << ",";
	std::cout << "]" << std::endl;
#endif

	return 0;
}
