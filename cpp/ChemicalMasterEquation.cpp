//g++ ChemicalMasterEquation.cpp -O3 -std=c++11 `pkg-config opencv --cflags --libs` -o CME
#include <iostream>
#include <random>
#include <opencv2/opencv.hpp>
#include <opencv2/plot.hpp>


void BrusselatorCME (std :: vector < double > & x, std :: vector < double > & y,
                     std :: vector < double > & t,
                     const double & A, const double & B,
                     const double & omega, double max_time = 30.,
                     std :: size_t seed = 123)
{
  double ti = t[0];
  double xi = x[0];
  double yi = y[0];

  std :: mt19937 mt(seed);
  std :: uniform_real_distribution < double > uniform(0., 1.);

  while (ti <= max_time)
  {
    std :: cout << "\rTime " << ti << std :: flush;

    std :: array < double, 4 > c{ {A * omega,
                                   A * omega + B * xi,
                                   A * omega + B * xi + xi,
                                   A * omega + B * xi + xi + xi * (xi - 1.) * yi / (omega * omega)
                                   }
                                 };

    const double rng1 = uniform(mt);
    const double rng2 = uniform(mt);
    const double tau = -std :: log(rng1) / c[3];
    const double uct = rng2 * c.back();

    ti += tau;

    if     (uct < c[0]) ++xi;
    else if(uct < c[1]) {--xi; ++yi;}
    else if(uct < c[2]) --xi;
    else if(uct < c[3]) {++xi; --yi;}

    x.push_back(xi);
    y.push_back(yi);
    t.push_back(ti);
  }

}

int main (int argc, char ** argv)
{
  const double A = 2.;
  const double B = 5.2;

  const double omega = 1000.;

  const double x0 = 1.6;
  const double y0 = 2.8;

  std :: vector < double > x = {x0};
  std :: vector < double > y = {y0};
  std :: vector < double > t = {0.};

  BrusselatorCME (x, y, t, A, B, omega, 30., 42);

  cv :: Mat plot_x;
  cv :: Mat plot_y;

  cv :: Ptr < cv :: plot :: Plot2d > plot = cv :: plot :: Plot2d :: create(cv :: Mat(x));
  plot->setPlotBackgroundColor( cv :: Scalar( 0, 0, 0 ) );
  plot->setPlotLineColor( cv :: Scalar( 255, 0, 0 ) );
  plot->setPlotAxisColor( cv :: Scalar( 255, 255, 255 ) );
  plot->setPlotGridColor( cv :: Scalar( 127, 127, 127 ) );
  plot->setPlotLineWidth(2);
  plot->setInvertOrientation(true);
  plot->setShowText(false);

  plot->render( plot_x );

  cv :: Mat plot_res;

  plot = cv :: plot :: Plot2d :: create(cv :: Mat(y));
  plot->setPlotAxisColor( cv :: Scalar( 255, 255, 255 ) );
  plot->setPlotLineColor( cv :: Scalar( 0, 0, 255 ) );
  plot->setPlotGridColor( cv :: Scalar( 127, 127, 127 ) );
  plot->setPlotLineWidth(2);
  plot->setInvertOrientation(true);
  plot->setShowText(false);
  plot->render( plot_y );

  plot_res = plot_x | plot_y;

  cv :: namedWindow("Brusselator CME", cv :: WINDOW_FULLSCREEN);
  cv :: moveWindow("Brusselator CME", 0, 10);
  cv :: imshow("Brusselator CME", plot_res);

  cv :: waitKey(0);

  return 0;
}
