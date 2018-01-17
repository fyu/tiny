#include "furry/common/io.h"

#include <iomanip>

#include <boost/io/ios_state.hpp>

namespace furry {

void BeginSession(const std::string &name,
                  int title_size,
                  int num_delimiters,
                  char delimiter) {
  int margin = (title_size - num_delimiters * 2 - name.size()) / 2;

  boost::io::ios_all_saver guard(std::cout);
  std::cout << std::setfill(delimiter) << std::setw(num_delimiters - 1)
            << delimiter
            << std::setfill(' ')
            << std::setw(margin + name.size()) << std::right << name
            << std::setw(title_size - 2 * num_delimiters -
                         margin - name.size() - 1)
            << ' '
            << std::setfill(delimiter) << std::setw(num_delimiters - 1)
            << delimiter << '\n';
}

void EndSession(const std::string &name,
                int title_size,
                int num_delimiters,
                char delimiter) {
  BeginSession(name, title_size, num_delimiters, delimiter);
}

} // furry
