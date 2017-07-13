#include "furry/common/cast.h"

#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>

namespace furry {

std::vector<int>
ParseIntListString(const std::string &str) {
  std::vector<std::string> split_vector;
  boost::algorithm::split(split_vector, str, boost::algorithm::is_any_of(", "),
                          boost::algorithm::token_compress_on);
  std::vector<int> result;
  for (auto &s : split_vector) {
    result.push_back(to<int, std::string>(s));
  }
  return result;
}

} // furry
