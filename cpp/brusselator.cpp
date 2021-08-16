//g++ brusselator.cpp -O3 -std=c++11 `pkg-config opencv --cflags --libs` -o brusselator -lpthread
#include <iostream>
#include <thread>
#include <opencv2/opencv.hpp>


void view (const std :: string & name, cv :: Mat & img, int32_t ms=1)
{
  cv :: Mat temp = img.clone();
  cv :: normalize(img, temp, 0, 255, cv :: NORM_MINMAX);
  temp.convertTo(temp, CV_8UC1 );

  cv :: applyColorMap(temp, temp, cv::COLORMAP_JET);
  cv :: imshow(name, temp);
  int c = cv :: waitKey(ms);
  c = (c != -1) ? c % 256 : c;

  if (c == 27)
  {
    cv :: destroyAllWindows();

    if (ms == 0)
      return;

    //std :: cout << std :: endl;
    std :: exit(0);
  }
}


void diffusion (cv :: Mat & U, cv :: Mat & V,
                const double & dt,
                const double & A, const double & B,
                const double & Du, const double & Dv,
                const int64_t & iteration)
{

  const std :: string name = "Turing Pattern";
  cv :: namedWindow(name, cv :: WINDOW_FULLSCREEN );


  for (int64_t t = 0; t < iteration; ++t)
  {
    std :: thread display = std :: thread(view, name, std :: ref(U), 1);

    cv :: Mat lap_u;
    cv :: Mat lap_v;

    // a wrap boundary condition should be more appropriated
    // but unfortunately OpenCV Laplacian doesn't support it :(
    cv :: Laplacian(U, lap_u, CV_64FC1, 1, 1, 0, cv :: BORDER_REFLECT);
    cv :: Laplacian(V, lap_v, CV_64FC1, 1, 1, 0, cv :: BORDER_REFLECT);

    cv :: Mat uuv = U.mul(U.mul(V));

    cv :: Mat Ut = dt * (Du * lap_u + (A - (B + 1.) * U + uuv)) + U;
    cv :: Mat Vt = dt * (Dv * lap_v + (B * U - uuv)) + V;

    U = Ut.clone();
    V = Vt.clone();

    //std :: cout << "\rTime: " << dt * t << std :: flush;
    cv :: setWindowTitle(name, name + " (Time: " + std :: to_string(dt * t) + ")");
    display.join();
  }
  //std :: cout << "\rTime: " << dt * iteration << std :: endl;
}


void usage (char ** argv)
{
  std :: cerr << "Usage: " << argv[0] << " [A <double>] [B <double>] [Dx <double>] [Dy <double>] [dt <double>]"
              << std :: endl
              << "Default parameters:" << std :: endl
              << "\tA = 4.5" << std :: endl
              << "\tB = 4.75" << std :: endl
              << "\tDx = 2.0" << std :: endl
              << "\tDy = 16.0" << std :: endl
              << "\tdt = 0.005" << std :: endl
              << std :: endl;
  std :: exit(1);
}




void parse_args (int32_t argc, char ** argv,
                 double & A, double & B, double & Dx, double & Dy,
                 double & dt)
{
  switch (argc)
  {
    default: usage(argv);

    case 1:
    break;

    case 2:
    {
      A = std :: stod(argv[1]);
    } break;
    case 3:
    {
      B = std :: stod(argv[2]);
      parse_args(2, argv, A, B, Dx, Dy, dt);
    } break;
    case 4:
    {
      Dx = std :: stod(argv[3]);
      parse_args(3, argv, A, B, Dx, Dy, dt);
    } break;
    case 5:
    {
      Dy = std :: stod(argv[4]);
      parse_args(4, argv, A, B, Dx, Dy, dt);
    } break;
    case 6:
    {
      dt = std :: stod(argv[5]);
      parse_args(5, argv, A, B, Dx, Dy, dt);
    } break;
  }
}


int main (int argc, char ** argv)
{
  const int64_t dim = 512;

  double A = 4.5;
  double B = 4.5;
  double Du = 2.;
  double Dv = 16.;
  double dt = .005;

  parse_args(argc, argv, A, B, Du, Dv, dt);

  cv :: Mat U(dim, dim, CV_64FC1);
  cv :: Mat V(dim, dim, CV_64FC1);

  cv :: randu(U, cv :: Scalar(0.), cv :: Scalar(1.));
  cv :: randu(V, cv :: Scalar(0.), cv :: Scalar(1.));

  U = A + .3 * U;
  V = B/A + .3 * V;

  view("Initial condition", U, 0);

  diffusion(U, V, dt, A, B, Du, Dv, 6000);

  return 0;
}
