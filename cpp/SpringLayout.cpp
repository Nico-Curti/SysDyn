#include <iostream>
#include <algorithm>
#include <unordered_set>
#include <iterator>
#include <numeric>
#include <cmath>
#include <random>
#include <array>
#include <memory>

static constexpr int DIM = 3;
static constexpr int SEED = 123;

std::mt19937 eng(SEED);
std::uniform_real_distribution<float> distr(0.1f, .2f);

struct node
{
  std::array<float, DIM> force, pos;
  node()
  {
    std::fill_n(this->force.begin(), DIM, 0.f);
    std::fill_n(this->pos.begin(), DIM, 0.f);
  };
  ~node() = default;
};

struct
{
  void operator()(node &n1,
                  node &n2,
                  const float &k,
                  const float &r)
  {
    std::array<float, DIM> delta;
    std::transform(n1.pos.begin(), n1.pos.end(),
                   n2.pos.begin(), delta.begin(),
                   [](const float &x1, const float &x2)
                   {
                    return x2 - x1;
                   });
    float distance = sqrt(std::accumulate(delta.begin(), delta.end(),
                                          0.f, [](const float &r, const float &d)
                                          {
                                            return r + d*d;
                                          }));
    // If the deltas are too small, use random values to keep things moving
    if(distance < .1f)
    {
      std::generate_n(delta.begin(), DIM, [](){return distr(eng);});
      distance = sqrt(std::accumulate(delta.begin(), delta.end(),
                                        0.f, [](const float &r, const float &d)
                                        {
                                          return r + d*d;
                                        }));
    }

    if(distance < r)
    {
      float force = (k*k)/(distance*distance);
      std::transform(n1.force.begin(), n1.force.end(),
                     delta.begin(), n1.force.begin(),
                     [&force](const float &f, const float &d)
                     {
                        return f - force*d;
                     });
      std::transform(n2.force.begin(), n2.force.end(),
                     delta.begin(), n2.force.begin(),
                     [&force](const float &f, const float &d)
                     {
                        return f + force*d;
                     });
    }
    return;
  }

} Coulomb;

struct
{
  void operator()(node &n1,
                  node &n2,
                  const float &k,
                  const float &r)
  {
    std::array<float, DIM> delta;
    std::transform(n1.pos.begin(), n1.pos.end(),
                   n2.pos.begin(), delta.begin(),
                   [](const float &x1, const float &x2)
                   {
                    return x2 - x1;
                   });
    float distance = std::sqrt(std::accumulate(delta.begin(), delta.end(),
                                               0.f, [](const float &r, const float &d)
                                               {
                                                 return r + d*d;
                                               }));
    // If the deltas are too small, use random values to keep things moving
    if(distance < .1f)
    {
      std::generate_n(delta.begin(), DIM, [](){return distr(eng);});
      distance = std::sqrt(std::accumulate(delta.begin(), delta.end(),
                                           0.f, [](const float &r, const float &d)
                                           {
                                             return r + d*d;
                                           }));
    }

    // Truncate distance so as to not have crazy springiness
    distance = std::min(distance, r);

    // Calculate Hooke force and update nodes
    float force = (distance*distance - k*k) / (distance*k);
    std::transform(n1.force.begin(), n1.force.end(),
                   delta.begin(), n1.force.begin(),
                   [&force](const float &f, const float &d)
                   {
                      return f + force*d;
                   });
    std::transform(n2.force.begin(), n2.force.end(),
                   delta.begin(), n2.force.begin(),
                   [&force](const float &f, const float &d)
                   {
                      return f - force*d;
                   });
    return;
  }

} Hooke;


struct
{
  auto operator()(std::pair<int, int> * edges,
                  const int &Nedges,
                  int iterations = 1000,
                  float force_strength = 5.f,
                  float damping = .01f,
                  float max_velocity = 2.f,
                  float max_distance = 50.f
                  )
  {
    std::unordered_set<int> n;
    for(int i = 0; i < Nedges; ++i)
    {
      n.insert(edges[i].first);
      n.insert(edges[i].second);
    }
    int i = 0, Nnodes = static_cast<int>(n.size());
    std::shared_ptr<node[]> nodes(new node[Nnodes]);
    while( i < iterations)
    {
      // Add in Coulomb-esque node-node repulsive forces
      for(int n1 = 0; n1 < Nnodes; ++n1)
        for(int n2 = 0; n2 < n1; ++n2)
          Coulomb(nodes[n1], nodes[n2], force_strength, max_distance);
      // And Hooke-esque edge spring forces
      for(int e = 0; e < Nedges; ++e)
        Hooke(nodes[edges[e].first], nodes[edges[e].second], force_strength, max_distance);
      // Move by resultant force
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for(int n = 0; n < Nnodes; ++n)
      {
        std::transform(nodes[n].force.begin(), nodes[n].force.end(),
                       nodes[n].pos.begin(), nodes[n].pos.begin(),
                       [&max_velocity, &damping](const float &f, const float &p)
                       {
                         return p + std::max(-max_velocity, std::min(damping*f, max_velocity));
                       });
        std::fill_n(nodes[n].force.begin(), DIM, 0.f);
      }
      ++i;
    }
    // Clean and return
    return nodes;
  }

} spring_layout;

int main(int argc, char **argv)
{
    /*======================================

        (0)-------(1)   Adjacency Matrix
        / \        /        0 1 1 1 0
       /   \      /         1 0 0 0 1
      (3)--(2)   /          1 0 0 1 0
        \       /           1 0 1 0 1
         \____(4)           0 1 0 1 0

    n_link = 6

    ======================================*/
    const int Nnodes = 5, Nedges = 12;
    std::array<std::pair<int, int>, Nedges> edges;
    edges[0]  = std::make_pair(0, 1);
    edges[1]  = std::make_pair(0, 2);
    edges[2]  = std::make_pair(0, 3);
    edges[3]  = std::make_pair(1, 4);
    edges[4]  = std::make_pair(3, 2);
    edges[5]  = std::make_pair(3, 4);
    edges[6]  = std::make_pair(1, 0);
    edges[7]  = std::make_pair(2, 0);
    edges[8]  = std::make_pair(3, 0);
    edges[9]  = std::make_pair(4, 1);
    edges[10] = std::make_pair(2, 3);
    edges[11] = std::make_pair(4, 3);

    auto nodes = spring_layout(edges.data(), Nedges);
#ifdef DEBUG
    std::cout << "pos = {" << std::endl;
    for(int i = 0; i < Nnodes; ++i)
    {
      std::cout << i << ":[";
      for(int j = 0; j < DIM-1; ++j) std::cout << nodes[i].pos[j] << ", ";
      std::cout << nodes[i].pos[DIM-1] << "]," << std::endl;
    }
    std::cout << "}\n";
#else
    std::cout << "Nodes\tx\ty\tz" << std::endl;
    for(int i = 0; i < Nnodes; ++i)
    {
      std::cout << i << "\t";
      std::copy_n(nodes[i].pos.begin(), DIM, std::ostream_iterator<float>(std::cout, "\t"));
      std::cout << std::endl;
    }
#endif

    return 0;
}