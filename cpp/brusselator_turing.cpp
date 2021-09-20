//g++ brusselator.cpp -O3 -std=c++11 `pkg-config opencv --cflags --libs` -o brusselator -lpthread
#include <iostream>
#include <thread>
#include <opencv2/opencv.hpp>

/**
* @brief OpenCV viewer
*
* @details This function is used for the visualization of
* an OpenCV image and it could be used for an asynchronous
* visualization of the result.
*
* @note The input image is normalized between its Min-Max
* and converted to uint8_t before the visualization.
* A Jet colormap (do not tell to prof. Giampieri that I have
* used a Jet colormap, please!) is used for the color remapping.
*
* @param name Window name
* @param img OpenCV Mat to plot
* @param ms Wait time in ms
*
*/
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

/**
* @brief Brusselator diffusion model
*
* @note The visualization of the 1st morphogen
* is performed asynchronously during the update
* computation using std :: thread.
*
* @param U OpenCV Mat of the 1st morphogen
* @param V OpenCV Mat of the 2nd morphogen
* @param dt Interval of time
* @param A Kinetic reaction constant
* @param B Kinetic reaction constant
* @param Du Diffusion coef of the 1st morphogen
* @param Dv Diffusion coef of the 2nd morphogen
* @param iteration Number of iterations to perform
*
*/
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
    cv :: Laplacian(U, lap_u, CV_64FC1, 1, 1, 0, cv :: BORDER_REFLECT_101);
    cv :: Laplacian(V, lap_v, CV_64FC1, 1, 1, 0, cv :: BORDER_REFLECT_101);

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


/**
* @brief Command line helper
*
* @details Utility function for the command line user.
*
* @param argv Array of command line arguments.
*
*/
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



/**
* @brief Command line parser
*
* @details Parse the command line arguments
* and (eventually) set default values of the
* required variables.
* If something goes wrong the helper function
* is called.
*
* @param argc Number of arguments in command line.
* @param argv Array of command line arguments.
* @param A Kinetic reaction constant
* @param B Kinetic reaction constant
* @param Dx Diffusion coef of the 1st morphogen
* @param Dy Diffusion coef of the 2nd morphogen
* @param dt Interval of time
*
*/
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
