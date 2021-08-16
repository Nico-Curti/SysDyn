#include <iostream>
#include <unordered_set>
#include <numeric>
#include <random>
#include <opencv2/opencv.hpp>

static const int32_t DIM = 2;
static const int32_t SEED = 42;

std :: mt19937 eng(SEED);
std :: uniform_real_distribution < float > distr(0.1f, .2f);

struct node
{
  std :: array < float, DIM > force;
  std :: array < float, DIM > pos;

  node ()
  {
    std :: fill_n(this->force.begin(), DIM, 0.f);
    std :: fill_n(this->pos.begin(), DIM, 0.f);
  };

  ~node() = default;
};

void view (const std :: string & name, const cv :: Mat & edges, node * nodes, const int32_t & Nnodes, int32_t ms=1)
{
  cv :: Mat pos(Nnodes, 2, CV_32FC1);
  for (int32_t n = 0; n < Nnodes; ++n)
    pos.at < cv :: Point2f >(n) = cv :: Point2f(nodes[n].pos[0], nodes[n].pos[1]);

  cv :: normalize(pos, pos, 6.f, 506.f, cv :: NORM_MINMAX);

  cv :: Mat canvas = cv :: Mat :: zeros(cv :: Size(512, 512), CV_8UC1);
  for (int32_t i = 0; i < Nnodes; ++i)
    cv :: circle(canvas, pos.at< cv :: Point2f >(i), 5, cv :: Scalar :: all(255), cv :: FILLED);

  for (int32_t i = 0; i < edges.rows; ++i)
  {
    cv :: Point2i edge = edges.at< cv :: Point2i >(i);
    cv :: Point2f start = pos.at < cv :: Point2f >(edge.x);
    cv :: Point2f end = pos.at < cv :: Point2f >(edge.y);
    cv :: line(canvas, start, end, cv :: Scalar :: all(128), 1);
  }

  cv :: imshow(name, canvas);
  int32_t c = cv :: waitKey(ms);
  c = (c != -1) ? c % 256 : c;

  if (c == 27)
  {
    cv :: destroyAllWindows();

    if (ms == 0)
      return;

    std :: exit(0);
  }
}




struct
{
  void operator()(node & n1,
                  node & n2,
                  const float & k,
                  const float & r)
  {
    std :: array < float, DIM > delta;

    std :: transform(n1.pos.begin(), n1.pos.end(),
                    n2.pos.begin(), delta.begin(),
                    [](const float & x1, const float & x2)
                    {
                     return x2 - x1;
                    });
    float distance = std :: accumulate(delta.begin(), delta.end(),
                                       0.f, [](const float & r, const float & d)
                                       {
                                         return r + d * d;
                                       });

    // If the deltas are too small, use random values to keep things moving
    if (distance < .1f)
    {
      std :: generate_n(delta.begin(), DIM, [](){return distr(eng);});
      distance = std::accumulate(delta.begin(), delta.end(),
                                 0.f, [](const float & r, const float & d)
                                 {
                                   return r + d * d;
                                 });
    }

    if (distance < r*r)
    {
      float force = (k * k) / distance;
      std :: transform(n1.force.begin(), n1.force.end(),
                       delta.begin(), n1.force.begin(),
                       [&](const float & f, const float & d)
                       {
                          return f - force * d;
                       });
      std :: transform(n2.force.begin(), n2.force.end(),
                       delta.begin(), n2.force.begin(),
                       [&](const float & f, const float & d)
                       {
                          return f + force*d;
                       });
    }
    return;
  }

} Coulomb;

struct
{
  void operator()(node & n1,
                  node & n2,
                  const float & k,
                  const float & r)
  {
    std :: array < float, DIM > delta;
    std :: transform(n1.pos.begin(), n1.pos.end(),
                     n2.pos.begin(), delta.begin(),
                     [](const float & x1, const float & x2)
                     {
                      return x2 - x1;
                     });
    float distance = std :: accumulate(delta.begin(), delta.end(),
                                       0.f, [](const float & r, const float & d)
                                       {
                                         return r + d * d;
                                       });

    // If the deltas are too small, use random values to keep things moving
    if (distance < .1f)
    {
      std :: generate_n(delta.begin(), DIM, [](){return distr(eng);});
      distance = std::accumulate(delta.begin(), delta.end(),
                                 0.f, [](const float & r, const float & d)
                                 {
                                   return r + d * d;
                                 });
    }

    // Truncate distance so as to not have crazy springiness
    distance = std :: min(distance, r);

    // Calculate Hooke force and update nodes
    float force = (distance - k * k) / (distance * distance * k);
    std :: transform(n1.force.begin(), n1.force.end(),
                     delta.begin(), n1.force.begin(),
                     [&](const float & f, const float & d)
                     {
                        return f + force * d;
                     });
    std :: transform(n2.force.begin(), n2.force.end(),
                     delta.begin(), n2.force.begin(),
                     [&](const float & f, const float & d)
                     {
                        return f - force * d;
                     });
    return;
  }

} Hooke;


struct
{
  auto operator()(const cv :: Mat & edges,
                  int iterations = 1000,
                  float force_strength = 5.f,
                  float damping = .01f,
                  float max_velocity = 2.f,
                  float max_distance = 50.f
                  )
  {
    std :: unordered_set < int32_t > n;

    for (int32_t i = 0; i < edges.rows; ++i)
    {
      cv :: Point2i edge = edges.at< cv :: Point2i >(i);
      n.insert(edge.x);
      n.insert(edge.y);
    }

    int32_t i = 0;
    int32_t Nnodes = static_cast < int32_t >(n.size());

    std :: unique_ptr < node[] > nodes(new node[Nnodes]);

    const std :: string name = "Spring Layout";

    while ( i < iterations )
    {

      // Add in Coulomb-esque node-node repulsive forces
      for (int32_t n1 = 0; n1 < Nnodes; ++n1)
        for (int32_t n2 = 0; n2 < n1; ++n2)
          Coulomb(nodes[n1], nodes[n2], force_strength, max_distance);

      // And Hooke-esque edge spring forces
      for (int32_t e = 0; e < edges.rows; ++e)
      {
        cv :: Point2i edge = edges.at< cv :: Point2i >(e);
        Hooke(nodes[edge.x], nodes[edge.y], force_strength, max_distance);
      }

      // Move by resultant force
      for (int32_t n = 0; n < Nnodes; ++n)
      {
        std :: transform(nodes[n].force.begin(), nodes[n].force.end(),
                         nodes[n].pos.begin(), nodes[n].pos.begin(),
                         [&](const float & f, const float & p)
                         {
                           return p + std :: max(-max_velocity, std :: min(damping*f, max_velocity));
                         });

        std :: fill_n(nodes[n].force.begin(), DIM, 0.f);
      }
      ++i;


      view (name, edges, nodes.get(), Nnodes, 1);
      cv :: setWindowTitle(name, name + " (Iter: " + std :: to_string(i) + ")");
    }

    // Clean and return
    cv :: Mat pos(Nnodes, 2, CV_32FC1);
    for (int32_t n = 0; n < Nnodes; ++n)
      pos.at < cv :: Point2f >(n) = cv :: Point2f(nodes[n].pos[0], nodes[n].pos[1]);

    return pos;
  }

} spring_layout;

int main (int argc, char ** argv)
{
  const int Nnodes = 50, Nedges = 100;

  cv :: Mat edges(Nedges, 2, CV_32SC1);
  cv :: randu(edges, cv :: Scalar(0.), cv :: Scalar(Nnodes));

  cv :: Mat pos = spring_layout(edges);

  return 0;
}
