#ifndef _FURRY_COMMON_MAT_H_
#define _FURRY_COMMON_MAT_H_

#include "opencv2/core/core.hpp"

namespace furry
{

std::ostream& operator << (std::ostream& os, const cv::Point2d& p);
std::ostream& operator << (std::ostream& os, const cv::Point3d& p);
template<typename T, int n>
std::ostream& operator << (std::ostream& os, const cv::Vec<T, n>& p);

////////////////////////////////////////////////////////////
// inline implementation
////////////////////////////////////////////////////////////

inline std::ostream& operator << (std::ostream& os, const cv::Point2d& p)
{
  os << p.x << ' ' << p.y;
  return os;
}

inline std::ostream& operator << (std::ostream& os, const cv::Point3d& p)
{
  os << p.x << ' ' << p.y << ' ' << p.z;
  return os;
}

template<typename DataType, int n>
inline std::ostream& operator << (std::ostream& os, const cv::Vec<DataType, n>& v){
  os << v[0];
  for (int i = 1; i < n; ++i)
  {
    os << ' ' << v[i];
  }
  return os;
}

} // furry

#endif // _FURRY_COMMON_MAT_H_

