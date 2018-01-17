#ifndef FURRY_COMMON_POINT_SET
#define FURRY_COMMON_POINT_SET

#include <vector>

#include <opencv2/core/core.hpp>
#include <Eigen/Dense>

// #undef CGAL_CFG_NO_CPP0X_ARRAY
// #include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
// #include <CGAL/Point_set_2.h>

#include "furry/common/point.h"
#include "furry/algorithm/kdtree.h"

namespace furry
{

class PointSet2
{
public:
  PointSet2();
  PointSet2(const std::vector<cv::Point2d>& points);
  template<typename InputIterator>
  PointSet2(InputIterator first, InputIterator last) {
    insert(first, last);
  }
  // void insert(const cv::Point2d& p);
  void insert(const std::vector<cv::Point2d>& points);

  template <typename InputIterator>
  void insert(InputIterator first, InputIterator last) {
    typedef typename InputIterator::value_type PointType;
    static_assert(PointTraits<PointType>::kNumDimensions == 2,
                  "PointSet2 only accept 2d points");
    typename PointTraits<PointType>::X x;
    typename PointTraits<PointType>::Y y;
    // std::vector<Point_2> kpoints;
    // for (auto it = first; it != last; ++it) {
    //   kpoints.push_back(Point_2(x(*it), y(*it)));
    // }
    // point_set_.insert(kpoints.begin(), kpoints.end());
    std::vector<Eigen::Vector2d> points;
    for (auto it = first; it != last; ++it) {
      points.push_back(Eigen::Vector2d(x(*it), y(*it)));
    }
    point_set_.Build(points);
  }

  void reset();
  std::vector<cv::Point2d> range_search(cv::Rect rect);
  Eigen::Vector2d NearestNeighbor(const Eigen::Vector2d &p);
private:
  /* typedef CGAL::Exact_predicates_inexact_constructions_kernel  K; */
  /* typedef CGAL::Point_set_2<K>::Vertex_handle  VertexHandle; */
  /* typedef K::Point_2  Point_2; */

  /* CGAL::Point_set_2<CGAL::Exact_predicates_inexact_constructions_kernel> */
  KdTree<2, double> point_set_;
}; // PointSet

} // furry


#endif
