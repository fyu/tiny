#ifndef FURRY_COMMON_TYPE_TRAITS_H_
#define FURRY_COMMON_TYPE_TRAITS_H_

#include <type_traits>

namespace furry {

template <typename T>
struct remove_rcv {
  typedef typename std::remove_cv<
    typename std::remove_reference<T>::type>::type type;
};

} // furry

#endif // FURRY_COMMON_TYPE_TRAITS_H_
