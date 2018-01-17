#ifndef _FURRY_COMMON_RAND_H_
#define _FURRY_COMMON_RAND_H_

#include <vector>
#include <cstdlib>
#include <cmath>

#ifndef M_PI
#define M_PI 3.1415926536
#endif

namespace furry
{

// Return a random int in [0, n)
inline int rand_int(int n)
{
  return std::rand() % n;
}

// return a random number uniform in [begin, end)
inline double rand(double begin, double end)
{
  return std::rand() / (static_cast<double>(RAND_MAX) + 1.) * (end - begin)
    + begin;
}

// Return a random number with standard normal distribution
// Box Muller transform
inline double randn()
{
  double u1 = 1 - rand(0, 1);
  double u2 = 1 - rand(0, 1);
  return std::sqrt(-2 * std::log(u1)) * std::cos(2 * M_PI * u2);
}

inline std::vector<size_t>
sample_unique(const size_t n, const size_t m) {
  // Return an array of random sampled unique number
  // n specify the range and m specify the number of random numbers
  // The membership is random but the order is not random.
  // If the random order is also desired, you can permutate it randomly
  std::vector<size_t> arr;
  if (m > n)
  {
    for (size_t i = 0; i < n; ++i)
    {
      arr.push_back(i);
    }

    return arr;
  }
  size_t rm = m; // The remaining number of numbers sill needed
  size_t rn = n; // The remaining number of numbers from 0 to n-1
  size_t cn = 0; // Record the corrent number
  while (arr.size() < m) {
    // For each number, the probabiliby it is selected is rm / rn
    if ((size_t)rand_int(rn) < rm) {
      arr.push_back(cn);
      --rm;
    }
    --rn;
    ++cn;
  }
  return arr;
}

template <typename T> std::vector<T>
sample_unique(const std::vector<T>& arr, const size_t m)
{
  if (m > arr.size())
    return arr;
  size_t rm = m;
  size_t rn = arr.size();
  size_t cn = 0;
  std::vector<T> samples;
  while (samples.size() < m)
  {
    if ((size_t)rand_int(rn) < rm)
    {
      samples.push_back(arr[cn]);
      --rm;
    }
    --rn;
    ++cn;
  }
  return samples;
}

template <typename T> std::vector<T>
SampleUnique(const std::vector<T>& arr, const size_t m) {
  return sample_unique(arr, m);
}

inline std::vector<size_t>
SampleUnique(const size_t n, const size_t m) {
  return sample_unique(n, m);
}

// void SampleUnique(int n, int m, std::vector<int> *samples);
void SampleUnique(int n, int m, int *samples);

void deal(int n, int k, int *a);
std::vector<int> deal(int n, int m);

} // furry

#endif // _FURRY_COMMON_RAND_H_
