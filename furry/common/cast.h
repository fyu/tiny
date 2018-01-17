#ifndef FURRY_COMMON_CAST
#define FURRY_COMMON_CAST

#include <string>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <functional>
#include <iterator>
#include <iomanip>
#include <type_traits>

#include "furry/common/type_traits.h"

namespace furry
{

template <typename T, typename S>
T saturate_cast(const S& s);

template <> inline
unsigned char saturate_cast<unsigned char, float>(const float &f) {
  if (f < 0) return 0;
  if (f > 255) return 255;
  return f;
}

template <typename Enabled, typename Target, typename... Source>
struct ToImpl;

template <typename Target, typename... Source> Target
To(Source&&... arg) {
  return ToImpl<std::true_type, Target, typename remove_rcv<Source>::type...>::cast(
      std::forward<Source>(arg)...);
}

template <typename Target>
struct ToImpl<std::true_type, Target, Target> {
  static Target cast(const Target &arg) {
    return arg;
  }
};

template <template <typename, typename> class TargetContainer,
          template <typename, typename> class SourceContainer,
          typename Target,
          typename Source,
          typename TargetAllocator,
          typename SourceAllocator>
struct ToImpl<std::true_type, TargetContainer<Target, TargetAllocator>,
              SourceContainer<Source, SourceAllocator>> {
  static TargetContainer<Target, TargetAllocator>
      cast(const SourceContainer<Source, SourceAllocator> &values) {
    TargetContainer<Target, TargetAllocator> results;
    // results.reserve(values.size()); // doesn't work for list
    auto result = std::back_inserter(results);
    for (auto &value : values) {
      result = To<Target>(value);
    }
    return results;
  }
};

template <>
struct ToImpl<std::true_type, int, double> {
  static int cast(double d) {
    return d + 0.5;
  }
};

template <>
struct ToImpl<std::true_type, int, float> {
  static int cast(float f) {
    return f + 0.5;
  }
};

template <typename ArithmeticType>
struct ToImpl<std::true_type, typename std::enable_if<
                std::is_arithmetic<ArithmeticType>::value, std::string>::type,
              ArithmeticType> {
  static std::string cast(ArithmeticType n) {
    std::stringstream ss;
    ss << n;
    return ss.str();
  }
};

// Convert int to fixed width string
template <>
struct ToImpl<std::true_type, std::string, int, int, char> {
  static std::string cast(int number, int width, char padding) {
    std::stringstream ss;
    ss << std::setw(width) << std::setfill(padding) << number;
    return ss.str();
  }
};

template <>
struct ToImpl<std::true_type, std::string, int, int> {
  static std::string cast(int number, int width) {
    return ToImpl<std::true_type, std::string, int, int, char>::cast(number, width, '0');
  }
};

template <>
struct ToImpl<std::true_type, std::string, const char*> {
  static std::string cast(const char *const str) {
    return std::string(str);
  }
};

template <>
struct ToImpl<std::true_type, float, double> {
  static float cast(double x) {
    return x;
  }
};

template <>
struct ToImpl<std::true_type, double, float> {
  static double cast(float x) {
    return x;
  }
};

template <>
struct ToImpl<std::true_type, unsigned char, float> {
  static unsigned char cast(float x) {
    return saturate_cast<unsigned char>(x);
  }
};

template <>
struct ToImpl<std::true_type, int, std::string> {
  static int cast(const std::string &str) {
    return std::atoi(str.c_str());
  }
};

template <typename Target, typename Source> Target
to(const Source& arg);

template <>
inline double to<double, std::string>(const std::string & s)
{
  return std::atof(s.c_str());
}

// template <>
// inline double to<double, string>(std::string s)
// {
//   return std::atof(s.c_str());
// }

template <>
inline int to<int, std::string>(const std::string & s)
{
  return std::atoi(s.c_str());
}

// template <>
// inline int to<int, string>(std::string s)
// {
//   return std::atoi(s.c_str());
// }

template <>
inline std::string to<std::string, int>(const int& n)
{
  std::stringstream ss;
  ss << n;
  return ss.str();
}

// template <template <typename> class Container, typename T, typename S>
// struct To<Container<T>, Container<S>> {
//   static std::function<Container<T> (Container<S>)> f =
//       [] (const Container<S>& source) -> Container<T> {
//     Container<T> result;
//     auto output = std::back_inserter<result>;
//     for (auto &s : source) {
//       *output++ = to<T, S>(s);
//     }
//     return result;
//   };
// };

std::vector<int>
ParseIntListString(const std::string &s);

} // furry

#endif // FURRY_COMMON_CAST
