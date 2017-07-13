#include "fmeta.hpp"
#include <iostream>
#include <typeinfo>

using namespace std;
using namespace furry;

int main()
{
  //  typename true_::apply<int>::type::apply<double>::type x;
  //typename APPLY(APPLY(true_, int), double) x;
  //typename apply<apply<true_, int>::type, double>::type x;
  //typename apply2<false_, int, double>::type x;
  //typename first_::apply<true_, double>::type x;
  //typename false_::apply<int, double>::type x;
  //typename if_::f<and_::f<true_, false_>::type, int, double>::type x;
  typename if_::f<not_::f<false_>::type, int, double>::type x;
  cout << typeid(x).name() << endl;
  return 0;
}
