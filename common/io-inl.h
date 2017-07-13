#ifndef FURRY_COMMON_IO_H_
#error "This file can only be included by furry/common/io.h"
#endif

#include <iterator>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/uio.h>
#include <unistd.h>

#include <Eigen/Dense>

namespace furry {

// template <typename InputIterator> int
// Write(const std::string &file_name,
//       InputIterator first,
//       InputIterator last,
//       void *buffer,
//       size_t num_bytes) {
//   static const size_t buffer_size = 2048;
//   bool with_buffer = true;

//   typedef typename InputIterator::value_type ElemType;
  
// }

// template <typename PointType, typename OutputIterator> inline int
// ReadPoints(const std::string& file_name,
//            OutputIterator result,
//            void *buffer,
//            size_t num_bytes) {
//   static const size_t buffer_size = 1024;
//   bool with_buffer = true;

//   static_assert(PointType::ColsAtCompileTime == 1,
//                 //Eigen::internal::traits<PointType>::RowsAtCompileTime == 1,
//                 "The point type should be a one dimentional column matrix.");
//   static_assert(PointType::RowsAtCompileTime != Eigen::Dynamic,
//                 "Number of rows should be fixed.");
//   typedef typename PointType::Scalar Scalar;

//   size_t dims = PointType::RowsAtCompileTime;
//   size_t point_size = dims * sizeof(Scalar);
//   num_bytes -= num_bytes % point_size;

//   if (num_bytes == 0) {
//     with_buffer = false;
//     num_bytes = buffer_size * dims * sizeof(Scalar);
//     buffer = malloc(num_bytes);
//   }

//   int file = open(file_name.c_str(), O_RDONLY);
//   if (file < 0) {
//     printf("Error: Can't open file %s\n", file_name.c_str());
//     if (!with_buffer)
//       free(buffer);
//     return -1;
//   }
//   size_t read_bytes = 1;
//   size_t offset = 0;
//   unsigned char *char_buffer;

//   while (read_bytes > 0) {
//     char_buffer  = static_cast<unsigned char*>(buffer);
//     read_bytes = read(file, char_buffer + offset, num_bytes - offset);
//     while (read_bytes >= point_size) {
//       *result++ = PointType(reinterpret_cast<Scalar*>(char_buffer));
//       char_buffer += point_size;
//       read_bytes -= point_size;
//     }
//     memcpy(buffer, char_buffer, read_bytes);
//     offset = read_bytes;
//   }

//   if (!with_buffer)
//     free(buffer);

//   return 0;
// }

// template <typename InputIterator> inline int
// WritePoints(const std::string &file_name,
//             InputIterator first,
//             InputIterator last,
//             void *buffer,
//             size_t num_bytes) {
//   static const size_t buffer_size = 1024;
//   bool with_buffer = true;

//   typedef typename InputIterator::value_type PointType;

//   static_assert(PointType::ColsAtCompileTime == 1,
//                 //Eigen::internal::traits<PointType>::RowsAtCompileTime == 1,
//                 "The point type should be a one dimentional column matrix.");
//   static_assert(PointType::RowsAtCompileTime != Eigen::Dynamic,
//                 "Number of rows should be fixed.");

//   typedef typename PointType::Scalar Scalar;

//   size_t dims = PointType::RowsAtCompileTime;
//   size_t point_size = dims * sizeof(Scalar);
//   //asert(point_size == sizeof(PointType));
//   num_bytes -= num_bytes % point_size;

//   if (num_bytes == 0) {
//     with_buffer = false;
//     num_bytes = buffer_size * dims * sizeof(Scalar);
//     buffer = malloc(num_bytes);
//   }

//   int file = open(file_name.c_str(), O_WRONLY | O_CREAT, S_IRUSR | S_IWUSR);
//   if (file < 0) {
//     // printf("Error: Can't open file %s\n", file_name.c_str());
//     if (!with_buffer)
//       free(buffer);
//     return -1;
//   }

//   // printf("before first loop %d\n", last - first);
//   while (first != last) {
//     size_t bytes_to_write = 0;
//     unsigned char *char_buffer = static_cast<unsigned char*>(buffer);
//     while ((bytes_to_write < num_bytes) && (first != last)) {
//       memcpy(char_buffer + bytes_to_write, first->data(), point_size);
//       ++first;
//       bytes_to_write += point_size;
//     }
//     // printf("copy to buffer %u\n", bytes_to_write);
//     size_t bytes_written;
//     do {
//       bytes_written = write(file, char_buffer, bytes_to_write);
//       bytes_to_write -= bytes_written;
//       char_buffer += bytes_written;
//       // printf("write buffer to file %u\n", bytes_written);
//     } while (bytes_to_write);
//     // printf("write one loop %d\n", last - first);
//   }

//   close(file);

//   return 0;
// }

} // furry
