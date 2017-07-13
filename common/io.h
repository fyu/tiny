#ifndef FURRY_COMMON_IO_H_
#define FURRY_COMMON_IO_H_

#include <fstream>
#include <iostream>
#include <string>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/archive_exception.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/list.hpp>
#include <boost/type_traits.hpp>
#include <boost/filesystem.hpp>

#include "furry/common/dir.h"

namespace furry {

enum IoStatus {
  IO_FAIL = 0,
  IO_SUCCEED = 1
};

// template <typename InputIterator> int
// Write(const std::string &file_name,
//       InputIterator first,
//       InputIterator last,
//       void *buffer,
//       size_t num_bytes);

template <typename InputType> int
ReadText(const std::string &file_name,
         InputType &input,
         bool quiet = false) {
  std::ifstream ifs(AddExtension(file_name, ".txt"), std::ios_base::in);
  try {
    boost::archive::text_iarchive ar(ifs);
    ar >> input;
  } catch (boost::archive::archive_exception &e) {
    if (!quiet)
      std::cerr << "Error in reading " << AddExtension(file_name, ".txt")
                << " " << e.what() << '\n';
    ifs.close();
    return IO_FAIL;
  }
  ifs.close();
  return IO_SUCCEED;
}

template <typename InputType> int
ReadBinary(const std::string &file_name,
           InputType &input,
           bool quiet = false) {
  std::ifstream ifs(AddExtension(file_name, ".bin"),
                    std::ios_base::in | std::ios_base::binary);
  try {
    boost::archive::binary_iarchive ar(ifs);
    ar >> input;
  } catch (boost::archive::archive_exception &e) {
    if (!quiet)
      std::cerr << "Error in reading " << AddExtension(file_name, ".bin")
                << " " << e.what() << '\n';
    ifs.close();
    return IO_FAIL;
  }
  ifs.close();
  return IO_SUCCEED;
}

template <typename OutputType> int
WriteText(const std::string &file_name,
         //const typename remove_cvr<OutputType>::type &output) {
         const OutputType &output) {
  // std::cout << "writing text\n";
  std::ofstream ofs(AddExtension(file_name, ".txt"),
                    std::ios_base::out | std::ios_base::trunc);
  // std::cout << "open file\n";
  try {
    boost::archive::text_oarchive ar(ofs);
    // std::cout << "init archive\n";
    ar << output;
  } catch (boost::archive::archive_exception &e) {
    std::cerr << "Error in writing " << AddExtension(file_name, ".txt")
              << " " << e.what() << "\n";
    ofs.close();
    return IO_FAIL;
  }
  ofs.close();
  return IO_SUCCEED;
  // std::cout << "finish writing\n";
}

template <typename OutputType> int
WriteBinary(const std::string &file_name,
         const OutputType &output) {
  std::ofstream ofs(AddExtension(file_name, ".bin"),
                    std::ios_base::out | std::ios_base::trunc |
                    std::ios_base::binary);
  try {
    boost::archive::binary_oarchive ar(ofs);
    ar << output;
  } catch (boost::archive::archive_exception &e) {
    std::cerr << "Error in writing " << AddExtension(file_name, ".txt")
              << " " << e.what() << '\n';
    ofs.close();
    return IO_FAIL;
  }
  ofs.close();
  return IO_SUCCEED;
}

void BeginSession(const std::string &name,
                  int title_size = 40,
                  int num_delimiters = 10,
                  char delimiter = '*');

void EndSession(const std::string &name = "End",
                int title_size = 40,
                int num_delimiters = 10,
                char delimiter = '*');

} // furry

namespace boost {
namespace serialization {

// template <class Archive, class Tp>
// void serialize(Archive &ar, cv::Point_<Tp> &p, const unsigned int version) {
//   ar & p.x & p.y;
// }

// template <class Archive, class Tp>
// void serialize(Archive &ar, std::vector<cv::Point_<Tp>> &points,
//                const unsigned int version) {
//   ar & points;
//
//}
// template <class Archive, class T>
// void serialize(Archive &ar, std::vector<T> &v, const unsigned int version) {
//   ar & v;
//
// }
} // serialization
} // boost

#include "furry/common/io-inl.h"

#endif // FURRY_COMMON_IO_H_
