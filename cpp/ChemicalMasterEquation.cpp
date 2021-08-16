//g++ ChemicalMasterEquation.cpp -O3 -std=c++11 `pkg-config opencv --cflags --libs` -o CME
#include <iostream>
#include <opencv2/opencv.hpp>
#include <opencv2/plot.hpp>

static constexpr float A  = 2.f;
static constexpr float B  = 5.2f;
static constexpr float dt = .01f;
static constexpr int ITER = 10;

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

void reaction (const float & x, const float & y,
               float * w, const float & noise)
{
  w[0] = A * noise;
  w[1] = B * x;
  w[2] = x;
  w[3] = x * (x - 1) * y / (noise * noise);
}

void reac_step (float & x, float & y,
                const float * comp, const float & act_step)
{
  if     (act_step < comp[0]) ++x;
  else if(act_step < comp[1]) {--x; ++y;}
  else if(act_step < comp[2]) --x;
  else if(act_step < comp[3]) {++x; --y;}
}

void CME (std :: vector < float > & x, std :: vector < float > & y,
          const float & omega, std :: size_t seed = 123)
{
  float trans = 0.f;
  float t_old = 0.f;
  float t     = 0.f;
  float tmp_x = x[0];
  float tmp_y = y[0];

  std :: array < float, 4 > w;
  std :: array < float, 4 > c;
  std :: srand(seed);

  while (t <= ITER + trans)
  {

    std :: cout << "\rTime " << t << std :: flush;

    reaction(tmp_x, tmp_y, w.data(), omega);

    c[0] = w[0];

    for (int32_t j = 1; j < 4; ++j)
      c[j] = c[j-1] + w[j];

    const float z1 = static_cast < float > (std :: rand()) / static_cast < float > (RAND_MAX);
    const float z2 = static_cast < float > (std :: rand()) / static_cast < float > (RAND_MAX);
    const float tau = - std :: log(z1) / c[3];
    const float uct = z2 * c[3];

    t += tau;

    reac_step(tmp_x, tmp_y, c.data(), uct);

    if (t > trans && t > t_old + dt)
    {
      x.push_back(tmp_x);
      y.push_back(tmp_y);
      t_old = t;
    }
  }
}

int main (int argc, char ** argv)
{
  const float omega   = 1000.f;
  const float x_start = 1.6f;
  const float y_start = 2.8f;

  std :: vector < float > x = {x_start};
  std :: vector < float > y = {y_start};

  CME (x, y, omega);

  cv :: Mat plot_x;
  cv :: Mat plot_y;

  cv :: Ptr < cv :: plot :: Plot2d > plot = cv :: plot :: Plot2d :: create(cv :: Mat(x));
  plot->setPlotBackgroundColor( cv :: Scalar( 0, 0, 0 ) );
  plot->setPlotLineColor( cv :: Scalar( 255, 0, 0 ) );
  plot->setPlotAxisColor( cv :: Scalar( 255, 255, 255 ) );
  plot->setPlotLineWidth(2);
  plot->setInvertOrientation(true);
  plot->setShowText(false);

  plot->render( plot_x );

  cv :: Mat plot_res;

  plot = cv :: plot :: Plot2d :: create(cv :: Mat(y));
  plot->setPlotAxisColor( cv :: Scalar( 0, 0, 0 ) );
  plot->setPlotLineColor( cv :: Scalar( 0, 0, 255 ) );
  plot->setPlotLineWidth(2);
  plot->setInvertOrientation(true);
  plot->setShowText(false);
  plot->render( plot_y );

  plot_res = plot_x | plot_y;

  cv :: namedWindow("Brusselator CME", cv :: WINDOW_FULLSCREEN);
  cv :: moveWindow("Brusselator CME", 0, 10);
  cv :: imshow("Brusselator CME", plot_res);

  return 0;
}
