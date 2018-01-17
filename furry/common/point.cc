#include "furry/common/point.h"

#include <cstdio>
#include <cstdint>

#include "furry/common/memory.h"
#include "furry/common/log.h"

namespace furry {

// double Distance2(const Point2 &p0, const Point2 &p1) {
//   return (p0.x() - p1.x()) * (p0.x() - p1.x()) +
//       (p0.y() - p1.y()) * (p0.y() - p1.y());
// }

// std::ostream& operator << (std::ostream& os, const Point2 &p) {
//   os << p.x() << ' ' << p.y();
//   return os;
// }

// std::istream& operator >> (std::istream& is, Point2 &p) {
//   is >> p.x() >> p.y();
//   return is;
// }

// std::ostream& operator << (std::ostream& os, const Point3 &p) {
//   os << p.x() << ' ' << p.y() << ' ' << p.z();
//   return os;
// }

// std::istream& operator >> (std::istream& is, Point3 &p) {
//   is >> p.x() >> p.y() >> p.z();
//   return is;
// }

// bool operator == (const Point3 &p0, const Point3 &p1) {
//   return p0.x() == p1.x() && p0.y() == p1.y() && p0.z() == p1.z();
// }

std::istream& operator >> (std::istream& is, Point3 &p) {
  is >> p.x() >> p.y() >> p.z();
  return is;
}

template <typename T>
int WritePoints(const T* const buffer, int dim, int num_points, FILE *file) {
  if (fwrite(&dim, sizeof(int), 1, file) != 1)
    return 0;
  int type_size = sizeof(T);
  if (fwrite(&type_size, sizeof(int), 1, file) != 1)
    return 0;
  if (fwrite(&num_points, sizeof(int), 1, file) != 1)
    return 0;
  if (fwrite(buffer, dim * sizeof(T), num_points, file) != num_points)
    return 0;
  return 1;
}

template <typename T>
int WritePoints(const std::vector<cv::Point_<T>> &points, FILE *file) {
  auto buffer = ScopedArray<T>(new T[2 * points.size()]);
  int index = -1;
  for (auto p : points) {
    buffer[++index] = p.x;
    buffer[++index] = p.y;
  }
  return WritePoints(buffer.Get(), 2, points.size(), file);
}

template <typename T, int Dim>
int WritePoints(const std::vector<Eigen::Matrix<T, Dim, 1>> &points,
                FILE *file) {
  auto buffer = ScopedArray<T>(new T[Dim * points.size()]);
  for (size_t i = 0; i < points.size(); ++i) {
    memcpy(buffer.Get() + i * Dim, points[i].data(), Dim * sizeof(T));
  }
  return WritePoints(buffer.Get(), Dim, points.size(), file);
}

template <typename T>
int ReadPoints(FILE *file, T **buffer, int *dim, int *type_size,
               int *num_points) {
  if (fread(dim, sizeof(int), 1, file) != 1)
    return 0;
  if (fread(type_size, sizeof(int), 1, file) != 1)
    return 0;
  if (fread(num_points, sizeof(int), 1, file) != 1)
    return 0;
  *buffer = new T[*dim * *num_points];
  int count = 0;
  if ((count = fread(*buffer, sizeof(T), *num_points * *dim, file)) != *dim * *num_points) {
    delete *buffer;
    return 0;
  }

  return 1;
}

template <typename T>
int ReadPoints(FILE *file, std::vector<cv::Point_<T>> *points) {
  T *buffer;
  int dim, type_size, num_points;
  if (!ReadPoints(file, &buffer, &dim, &type_size, &num_points)) {
    F_ELOG("Failed to read buffer");
    return 0;
  }
  if (type_size != sizeof(T) || dim != 2)
    return 0;
  points->reserve(num_points + points->size());
  cv::Point_<T> point;
  int buffer_size = num_points * dim;
  for (int i = 0; i < buffer_size; i += 2) {
    point.x = buffer[i];
    point.y = buffer[i + 1];
    points->push_back(point);
  }
  delete [] buffer;
  return 1;
}

template <typename T, int Dim>
int ReadPoints(FILE *file, std::vector<Eigen::Matrix<T, Dim, 1>> *points) {
  T *buffer;
  int dim, type_size, num_points;
  if (!ReadPoints(file, &buffer, &dim, &type_size, &num_points))
    return 0;
  if (type_size != sizeof(T) || dim != Dim)
    return 0;
  points->reserve(num_points + points->size());
  Eigen::Matrix<T, Dim, 1> p;
  int buffer_size = num_points * dim;
  for (int i = 0; i < buffer_size; i += dim) {
    memcpy(p.data(), buffer + i, sizeof(T) * Dim);
    points->push_back(p);
  }
  delete [] buffer;
  return 1;
}

int WritePoints(const std::vector<CvPoint2> &points, FILE *file) {
  return WritePoints<Scalar2d>(points, file);
}

int ReadPoints(FILE *file, std::vector<CvPoint2> *points) {
  return ReadPoints<Scalar2d>(file, points);
}

int WritePoints(const std::vector<Vector2> &points, FILE *file) {
  return WritePoints<Scalar2d>(points, file);
}

int ReadPoints(FILE *file, std::vector<Vector2> *points) {
  return ReadPoints<Scalar2d>(file, points);
}

int WritePoints(const std::vector<Vector3> &points, FILE *file) {
  return WritePoints<Scalar3d>(points, file);
}

int ReadPoints(FILE *file, std::vector<Vector3> *points) {
  return ReadPoints<Scalar3d>(file, points);
}

int WritePoints(const std::vector<cv::KeyPoint> &points, FILE *file) {
  static const int point_size_in_byte = 28;
  int num_points = points.size();
  int num_bytes = point_size_in_byte * num_points;
  auto buffer =
      ScopedArray<uint8_t>(new uint8_t[num_bytes]);
  uint8_t *ptr = buffer.Get();
  for (int i = 0; i < num_points; ++i) {
    *((float*)ptr) = points[i].pt.x;
    *(((float*)ptr) + 1) = points[i].pt.y;
    *(((float*)ptr) + 2) = points[i].size;
    *(((float*)ptr) + 3) = points[i].angle;
    *(((float*)ptr) + 4) = points[i].response;
    *(((int*)ptr) + 5) = points[i].octave;
    *(((int*)ptr) + 6) = points[i].class_id;
    ptr += point_size_in_byte;
  }
  if (fwrite(&num_points, sizeof(int), 1, file) != 1) {
    return 1;
  }
  if (fwrite(buffer.Get(), 1, num_bytes, file) != num_bytes) {
    return 1;
  }
  return 0;
}

int ReadPoints(FILE *file, std::vector<cv::KeyPoint> *points) {
  static const int point_size_in_byte = 28;
  int num_points;
  if (fread(&num_points, sizeof(int), 1, file) != 1) {
    return 0;
  }
  int num_bytes = num_points * point_size_in_byte;
  auto buffer = ScopedPtr<uint8_t>(new uint8_t[num_bytes]);
  if (fread(buffer.Get(), 1, num_bytes, file) != num_bytes) {
    return 0;
  }
  uint8_t *ptr = buffer.Get();
  cv::KeyPoint point;
  for (int i = 0; i < num_points; ++i) {
    point.pt.x     = *((float*)ptr);
    point.pt.y     = *(((float*)ptr) + 1);
    point.size     = *(((float*)ptr) + 2);
    point.angle    = *(((float*)ptr) + 3);
    point.response = *(((float*)ptr) + 4);
    point.octave   = *(((int*)ptr) + 5);
    point.class_id = *(((int*)ptr) + 6);
    points->push_back(point);
    ptr += point_size_in_byte;
  }
  return 1;
}

int WritePoints(const std::string &filename,
                const std::vector<Point2> &points) {
  // return WritePoints(filename, points.begin(), points.end());
  auto file = ScopedPtr<FILE>(fopen(filename.c_str(), "w"));
  return WritePoints(points, file.Get());
}

int ReadPoints(const std::string &filename,
               std::vector<Point2> *points) {
  // return ReadPoints<Point2>(filename, std::back_inserter(*points));
  auto file = ScopedPtr<FILE>(fopen(filename.c_str(), "r"));
  return ReadPoints(file.Get(), points);
}

int WritePoints(const std::string &filename,
                const std::vector<Vector2> &points) {
  // return WritePoints(filename, points.begin(), points.end());
  auto file = ScopedPtr<FILE>(fopen(filename.c_str(), "w"));
  return WritePoints(points, file.Get());
}

int ReadPoints(const std::string &filename,
               std::vector<Vector2> *points) {
  // return ReadPoints<Vector2>(filename, std::back_inserter(*points));
  auto file = ScopedPtr<FILE>(fopen(filename.c_str(), "r"));
  return ReadPoints(file.Get(), points);
}

int WritePoints(const std::string &filename,
                const std::vector<Vector3> &points) {
  // return WritePoints(filename, points.begin(), points.end());
  auto file = ScopedPtr<FILE>(fopen(filename.c_str(), "w"));
  return WritePoints(points, file.Get());
}

int ReadPoints(const std::string &filename,
               std::vector<Vector3> *points) {
  // return ReadPoints<Vector3>(filename, std::back_inserter(*points));
  auto file = ScopedPtr<FILE>(fopen(filename.c_str(), "r"));
  return ReadPoints(file.Get(), points);
}

int ReadPoints(const std::string &filename,
               std::vector<cv::KeyPoint> *points) {
  auto file = ScopedPtr<FILE>(fopen(filename.c_str(), "r"));
  return ReadPoints(file.Get(), points);
}

int WritePoints(const std::string &filename,
                const std::vector<cv::KeyPoint> &points) {
  auto file = ScopedPtr<FILE>(fopen(filename.c_str(), "w"));
  return WritePoints(points, file.Get());
}



} // furry
