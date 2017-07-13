#include "meta.hpp"
#include <iostream>
#include <typeinfo>
#include <boost/mpl/if.hpp>
#include <boost/mpl/logical.hpp>

using namespace furry;
using namespace std;
using namespace boost;

//  template <typename T1, typename T2>
//  using t = typename to_bool_<false>::type<T1, T2>;

int main()
{

  typename if_<true_, int, double>::type x;
  //  typename mpl::if_<mpl::and_<mpl::true_, mpl::true_>, int, double>::type x;
  cout << typeid(x).name() << endl;
  return 0;
}
