#ifndef FURRY_COMMON_POINT_H_
#error "This file can only be included by furry/common/io.h"
#endif

#include <iterator>
#include <cmath>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/uio.h>
#include <unistd.h>
#include <iostream>

#include <Eigen/Dense>

#include "furry/common/log.h"

namespace furry {

template <typename PointType, typename OutputIterator> inline int
ReadPoints(const std::string& file_name,
           OutputIterator result,
           void *buffer,
           size_t num_bytes) {
  F_LOG(std::cout, "Warning: this function is deprecated");
  static const size_t buffer_size = 1024;
  bool with_buffer = true;

  // static_assert(PointType::ColsAtCompileTime == 1,
  //               //Eigen::internal::traits<PointType>::RowsAtCompileTime == 1,
  //               "The point type should be a one dimentional column matrix.");
  // static_assert(PointType::RowsAtCompileTime != Eigen::Dynamic,
  //               "Number of rows should be fixed.");
  typedef typename PointTraits<PointType>::ValueType Scalar;

  size_t dims = PointTraits<PointType>::kNumDimensions;
  size_t point_size = dims * sizeof(Scalar);
  num_bytes -= num_bytes % point_size;

  typename PointTraits<PointType>::SetData set_data;

  if (num_bytes == 0) {
    with_buffer = false;
    num_bytes = buffer_size * dims * sizeof(Scalar);
    buffer = malloc(num_bytes);
  }

  int file = open(file_name.c_str(), O_RDONLY);
  if (file < 0) {
    printf("Error: Can't open file %s\n", file_name.c_str());
    if (!with_buffer)
      free(buffer);
    return 0;
  }
  ssize_t read_bytes = 1;
  size_t offset = 0;
  unsigned char *char_buffer;

  while (true) {
    char_buffer  = static_cast<unsigned char*>(buffer);
    read_bytes = read(file, char_buffer + offset, num_bytes - offset);
    if (read_bytes == 0)
      break;
    if (read_bytes < 0) {
      printf("Error in reading file %s\n", file_name.c_str());
      break;
    }
    while (read_bytes >= (ssize_t)point_size) {
      //*result++ = PointType(reinterpret_cast<Scalar*>(char_buffer));
      PointType p;
      set_data(p, reinterpret_cast<Scalar*>(char_buffer));
      *result++ = p;
      char_buffer += point_size;
      read_bytes -= point_size;
    }
    memcpy(buffer, char_buffer, read_bytes);
    offset = read_bytes;
  }

  if (!with_buffer)
    free(buffer);

  return 1;
}

template <typename InputIterator> inline int
WritePoints(const std::string &file_name,
            InputIterator first,
            InputIterator last,
            void *buffer,
            size_t num_bytes) {
  F_LOG(std::cout, "Warning: this function is deprecated");
  static const size_t buffer_size = 1024;
  bool with_buffer = true;

  typedef typename InputIterator::value_type PointType;

  // static_assert(PointType::ColsAtCompileTime == 1,
  //               //Eigen::internal::traits<PointType>::RowsAtCompileTime == 1,
  //               "The point type should be a one dimentional column matrix.");
  // static_assert(PointType::RowsAtCompileTime != Eigen::Dynamic,
  //               "Number of rows should be fixed.");

  typedef typename PointTraits<PointType>::ValueType Scalar;
  typename PointTraits<PointType>::Data data;

  size_t dims = PointTraits<PointType>::kNumDimensions;
  size_t point_size = dims * sizeof(Scalar);
  //asert(point_size == sizeof(PointType));
  num_bytes -= num_bytes % point_size;

  if (num_bytes == 0) {
    with_buffer = false;
    num_bytes = buffer_size * dims * sizeof(Scalar);
    buffer = malloc(num_bytes);
  }

  int file = open(file_name.c_str(),
                  O_WRONLY | O_CREAT | O_TRUNC,
                  S_IRUSR | S_IWUSR);
  if (file < 0) {
    // printf("Error: Can't open file %s\n", file_name.c_str());
    if (!with_buffer)
      free(buffer);
    return 0;
  }

  // printf("before first loop %d\n", last - first);
  int num_points = 0;
  while (first != last) {
    size_t bytes_to_write = 0;
    unsigned char *char_buffer = static_cast<unsigned char*>(buffer);
    while ((bytes_to_write < num_bytes) && (first != last)) {
      memcpy(char_buffer + bytes_to_write, data(&(*first)), point_size);
      ++first;
      ++num_points;
      bytes_to_write += point_size;
    }
    // printf("copy to buffer %u\n", bytes_to_write);
    size_t bytes_written;
    do {
      bytes_written = write(file, char_buffer, bytes_to_write);
      bytes_to_write -= bytes_written;
      char_buffer += bytes_written;
      // printf("write buffer to file %u\n", bytes_written);
    } while (bytes_to_write);
    // printf("write one loop %d\n", last - first);
  }

  close(file);

  return num_points;
}

template <typename PointType0, typename PointType1> inline
typename std::enable_if<PointTraits<PointType0>::kNumDimensions ==
                        PointTraits<PointType1>::kNumDimensions,
                        typename PointTraits<PointType0>::ValueType>::type
EuclideanDistance(const PointType0 &p0, const PointType1 &p1) {
  return sqrt(EuclideanDistance2(p0, p1));
}

template <typename PointType0, typename PointType1> inline
typename std::enable_if<PointTraits<PointType0>::kNumDimensions ==
                        PointTraits<PointType1>::kNumDimensions,
                        typename PointTraits<PointType0>::ValueType>::type
EuclideanDistance2(const PointType0 &p0, const PointType1 &p1) {
  typename PointTraits<PointType0>::ValueType result = 0;
  typename PointTraits<PointType0>::Get get0;
  typename PointTraits<PointType1>::Get get1;
  for (size_t i = 0; i < PointTraits<PointType0>::kNumDimensions; ++i) {
    auto d = get0(p0, i) - get1(p1, i);
    result += d * d;
  }
  return result;
}

template <typename InputIterator, typename PointType> inline InputIterator
NearestNeighbor(InputIterator first, InputIterator last, const PointType &p) {
  if (first == last)
    return last;
  InputIterator nn = first;
  auto best_distance = furry::EuclideanDistance2(*first++, p);
  while (first != last) {
    auto distance = furry::EuclideanDistance2(*first, p);
    if (distance < best_distance) {
      best_distance = distance;
      nn = first;
    }
    ++first;
  }
  return nn;
}

} // furry
