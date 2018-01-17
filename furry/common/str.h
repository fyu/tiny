#ifndef FURRY_COMMON_STRING_H_
#define FURRY_COMMON_STRING_H_

#include <string>
#include <vector>
#include <sstream>
#include <cstdarg>

#include "furry/common/cast.h"

namespace furry
{

std::vector<std::string>
Split(const std::string & str, const std::string & delimiters);

std::string
tolower(std::string str);

template <typename ForwardIterator> std::string
Join(ForwardIterator first,
     ForwardIterator last,
     const std::string & delim)
{
  std::stringstream ss;
  if (first != last) ss << *first++;
  while (first != last)
  {
    ss << delim << *first++;
  }
  return ss.str();
}

namespace internal {

inline std::string Join(std::stringstream &ss, const std::string &delim) {
  return ss.str();
}

template <typename T>
std::string Join(std::stringstream &ss, const std::string &delim,
                 T&& value) {
  ss << value;
  return ss.str();
}

template <typename T1, typename T2, typename ... Args>
std::string Join(std::stringstream &ss, const std::string &delim,
                 T1&& v1, T2&& v2, Args&&... args) {
  ss << v1 << delim;
  return Join(ss, delim, v2, args...);
}

} // internal

template <typename ... Args>
std::string Join(const std::string & delim, Args&& ... args) {
  std::stringstream ss;
  return internal::Join(ss, delim, args...);
}

void StringAppendV(std::string* dst, const char* format, va_list ap);

std::string StringPrintf(const char* format, ...);

}

#endif // FURRY_COMMON_STRING
