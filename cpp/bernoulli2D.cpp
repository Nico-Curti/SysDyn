//g++ bernoulli.cpp -O3 -std=c++11 `pkg-config opencv --cflags --libs` -o bernoulli -lpthread
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


void bernoulli (cv :: Mat & G,
                const int64_t & iteration)
{

  const std :: string name = "Bernoulli Pattern";
  cv :: namedWindow(name, cv :: WINDOW_FULLSCREEN );


  for (int64_t t = 0; t < iteration; ++t)
  {
    std :: thread display = std :: thread(view, name, std :: ref(G), 10);

    G = (2 * G);
    cv :: subtract(G, 1., G, (G > 1.));

    //std :: cout << "\rIter: " << t << std :: flush;
    cv :: setWindowTitle(name, name + " (Iter: " + std :: to_string(t) + ")");
    display.join();
  }
  //std :: cout << "\rIter: " << iteration << std :: endl;
}



void usage (char ** argv)
{
  std :: cerr << "Usage: " << argv[0] << " [mx <double>] [my <double>] [sx <double>] [sy <double>]"
              << std :: endl
              << "Default parameters:" << std :: endl
              << "\tmx = 0.5" << std :: endl
              << "\tmy = 0.5" << std :: endl
              << "\tsx = 0.3" << std :: endl
              << "\tsy = 0.3" << std :: endl
              << std :: endl;
  std :: exit(1);
}


void parse_args (int32_t argc, char ** argv,
                 double & mx, double & my, double & sx, double & sy)
{
  switch (argc)
  {
    default: usage(argv);

    case 1:
    break;

    case 2:
    {
      mx = std :: stod(argv[1]);
    } break;
    case 3:
    {
      my = std :: stod(argv[2]);
      parse_args(2, argv, mx, my, sx, sy);
    } break;
    case 4:
    {
      sx = std :: stod(argv[3]);
      parse_args(2, argv, mx, my, sx, sy);
    } break;
    case 5:
    {
      sy = std :: stod(argv[4]);
      parse_args(2, argv, mx, my, sx, sy);
    } break;
  }
}


int main (int argc, char ** argv)
{
  const int64_t dim = 512;

  double mean_x = 0.5;
  double mean_y = 0.5;
  double std_x = 0.31;
  double std_y = 0.31;

  parse_args(argc, argv, mean_x, mean_y, std_x, std_y);

  cv :: Mat G(dim, dim, CV_64FC1);

  for (int32_t i = 0; i < dim; ++i)
    for (int32_t j = 0; j < dim; ++j)
    {
      const double x = static_cast < double >(i) / dim;
      const double y = static_cast < double >(j) / dim;
      const double gx = .5 * (x - mean_x)*(x - mean_x) / (std_x * std_x);
      const double gy = .5 * (y - mean_y)*(y - mean_y) / (std_y * std_y);
      G.at < double >(i, j) = std :: exp(- (gx + gy));
    }

  view("Initial condition", G, 0);

  bernoulli(G, 100);

  return 0;
}
