#ifndef FURRY_COMMON_FACTORY_H_
#define FURRY_COMMON_FACTORY_H_

#include <utility>
#include <functional>
#include <vector>

#include "furry/common/type_traits.h"

namespace furry {

template <typename Enabled, typename Product, typename... Args>
struct MakeImpl;

template <typename Product, typename... Args> Product
Make(Args&&... args) {
  return MakeImpl<std::true_type,
                  Product,
                  typename furry::remove_rcv<Args>::type...>::Do(
      std::forward<Args>(args)...);
}

template <typename T0, typename T1, typename T2, typename T3>
struct MakeImpl<std::true_type, std::pair<T0, T1>, std::pair<T2, T3>> {
  static std::pair<T0, T1> Do(const std::pair<T2, T3> &p) {
    return std::make_pair(Make<T0>(p.first), Make<T1>(p.second));
  }
};

template <typename Product, typename Arg>
struct MakeImpl<std::true_type, std::vector<Product>, std::vector<Arg>> {
  static std::vector<Product> Do(const std::vector<Arg> &args) {
    std::vector<Product> products;
    products.reserve(args.size());
    for (auto &arg : args) {
      products.push_back(Make<Product>(arg));
    }
    return products;
  }
};

} // furry

#endif
