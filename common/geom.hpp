#ifndef _FURRY_COMMON_GEOM_H_
#define _FURRY_COMMON_GEOM_H_

#include "mat.hpp"

namespace furry
{

template <typename T, int s, int t>
cv::Vec<T, t> project(const cv::Matx<T, t, s>& m, const cv::Vec<T, s>& v)
{
  return m * v;
}

// project homogeneous coordinates
template <typename T, int s, int t>
cv::Vec<T, t> project(const cv::Matx<T, t+1, s+1>& m, const cv::Vec<T, s>& v)
{
  cv::Mat hv(v);
  hv.push_back(1.);
  cv::Mat r = cv::Mat(m) * cv::Mat(hv);
  for (int i = 0; i < t; ++i)
  {
    r.at<T>(i, 0) /= r.at<T>(t, 0);
  }
  r.pop_back(1);
  return r;
}


} // furry


#endif // _FURRY_COMMON_GEOM_H_
