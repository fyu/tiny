#include "furry/common/point_set.h"
//#include <CGAL/assertions.h>

namespace furry
{

PointSet2::PointSet2()
{
}

PointSet2::PointSet2(const std::vector<cv::Point2d>& points)
{
  insert(points.begin(), points.end());
}

// void
// PointSet2::insert(const cv::Point2d& p)
// {
//   point_set_.insert(Eigen::Vector2d(p.x, p.y));
// }

void
PointSet2::insert(const std::vector<cv::Point2d>& points)
{
  insert(points.begin(), points.end());
}

void
PointSet2::reset()
{
  //  point_set_ = CGAL::Point_set_2<CGAL::Exact_predicates_inexact_constructions_kernel>();
  point_set_.Destroy();
}

std::vector<cv::Point2d>
PointSet2::range_search(cv::Rect rect)
{
  // std::vector<VertexHandle> handles;
  // Point_2 a(rect.x, rect.y + rect.height);
  // Point_2 b(rect.x, rect.y);
  // Point_2 c(rect.x + rect.width, rect.y);
  // Point_2 d(rect.x + rect.width, rect.y + rect.height);
  // point_set_.range_search(a, b, c, d, std::back_inserter(handles));
  // std::vector<cv::Point2d> points;
  // points.reserve(handles.size());
  // for (auto it = handles.begin(); it != handles.end(); ++it)
  // {
  //   auto p = (*it)->point();
  //   points.push_back(cv::Point2d(p.x(), p.y()));
  // }
  std::vector<Eigen::Vector2d> points;
  Eigen::Vector2d lower(rect.x, rect.y);
  Eigen::Vector2d upper(rect.x + rect.width, rect.y + rect.height);
  point_set_.RangeSearch(lower, upper, std::back_inserter(points));
  std::vector<cv::Point2d> cvpoints;
  for (auto it = points.begin(); it != points.end(); ++it) {
    cvpoints.push_back(cv::Point2d(it->x(), it->y()));
  }
  return cvpoints;
}

Eigen::Vector2d PointSet2::NearestNeighbor(const Eigen::Vector2d &p_) {
  // Point_2 p(p.x(), p.y());
  // auto h = point_set_.nearest_neighbor(p);
  // if (h != NULL) {
  //   auto nn = h->point();
  //   return Eigen::Vector2d(nn.x(), nn.y());
  // } else {
  //   return Eigen::Vector2d(-1, -1);
  // }
  return point_set_.NN(p_);
}

} // furry
