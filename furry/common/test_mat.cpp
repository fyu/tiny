#include "mat.hpp"
#include <iostream>

using namespace furry;
using namespace std;

int main()
{
  //cout << get_dims<1, 2, 3, 9, 10>::value << endl;
  // Vec<char, 5> d[10];
  // cout << sizeof(d) << endl;
  // cout << DataTraits<Vec<float, 5>>::type_size << " "
       // << DataTraits<Vec<char, 5>>::num_channels << endl;
  Vec<double, 3> v1;
  Vec<double, 3> v2;
  for (int i = 0; i < 3; ++i)
  {
    v1[i] = 1.5;
    v2[i] = 1.2;
  }
  v1 = v1 * 2;
  Vec<double, 3> v3 = v1 + v2;
  for (int i = 0; i < 3; ++i)
    cout << v3[i] << " ";
  cout << endl;
  return 0;
}
