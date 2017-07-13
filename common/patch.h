#ifndef FURRY_common_PATCH
#define FURRY_common_PATCH

#include <opencv2/opencv.hpp>

#include "furry/algorithm/kdtree.h"
#include "furry/common/cv.h"

namespace furry
{

// template <> inline cv::Rect
// to<cv::Rect>(double x,
//              double y,
//              double width,
//              double height)
// {
//   return cv::Rect(x - width / 2.0 + 0.5,
//                   y - height / 2.0 + 0.5,
//                   width + 0.5,
//                   height + 0.5);
// }

// template <> inline cv::Rect
// to<cv::Rect>(double x,
//              double y,
//              int width,
//              int height)
// {
//   return cv::Rect(x - width / 2.0 + 0.5,
//                   y - height / 2.0 + 0.5,
//                   width,
//                   height);
// }

// template <> inline cv::Rect
// to<cv::Rect>(int x,
//              int y,
//              int width,
//              int height)
// {
//   return cv::Rect(x - width / 2.0 + 0.5,
//                   y - height / 2.0 + 0.5,
//                   width,
//                   height);
// }

cv::Rect
to_rect(double x, double y, double width, double height);

cv::Rect
to_rect(const cv::Point2d& p, double s);

cv::Mat
make_patch(const cv::Mat& mat, double x, double y, double width, double height);

class PatchScoreFunc
{
public:
  cv::Scalar operator () (const cv::Mat&, const cv::Mat&) const;
  virtual bool better(double, double) const = 0;
  virtual double worst() const = 0;
  virtual ~PatchScoreFunc();

protected:
  virtual cv::Scalar compute(const cv::Mat&, const cv::Mat&) const = 0;
};

class PatchSsd : public PatchScoreFunc
{
public:
  virtual bool better(double, double) const;
  virtual double worst() const;
protected:
  virtual cv::Scalar compute(const cv::Mat&, const cv::Mat&) const;
};

class PatchNcc : public PatchScoreFunc
{
public:
  virtual bool better(double, double) const;
  virtual double worst() const;
protected:
  virtual cv::Scalar compute(const cv::Mat&, const cv::Mat&) const;
};

class PatchNccm : public PatchScoreFunc
{
public:
  virtual bool better(double, double) const;
  virtual double worst() const;
protected:
  virtual cv::Scalar compute(const cv::Mat&, const cv::Mat&) const;
};

cv::Point2d
adjust_patch_match(const cv::Mat& image0,
                   const cv::Mat& image1,
                   const cv::Point2d& p0,
                   const cv::Point2d& p1,
                   int patch_size,
                   int search_range,
                   PatchScoreFunc* psf,
                   double* final_score = NULL);

// match from set of patches to another set of patches in the search_range
// match quality is judged by psf
// result is a mapping from indices of centers0 to centers1
// If no match, -1 is assigned to the second of the pair
// std::vector<cv::Point2d>
// MatchProjectedFeatures(const cv::Mat& image0,
//                        const cv::Mat& image1,
//                        const std::vector<cv::Point2d> &points,
//                        const std::vector<cv::Point2d> &projected_points,
//                        const std::vector<cv::Point2d> &features,
//                        int patch_size,
//                        int search_range,
//                        PatchScoreFunc* psf,
//                        std::vector<double>* scores = NULL);

template <typename T> std::vector<cv::Point_<T>>
MatchProjectedFeatures(const cv::Mat& image0,
                       const cv::Mat& image1,
                       const std::vector<cv::Point_<T>> &points,
                       const std::vector<cv::Point_<T>> &projected_points,
                       const std::vector<cv::Point_<T>> &features,
                       int patch_size,
                       int search_range,
                       PatchScoreFunc* psf,
                       std::vector<double>* scores = NULL) {
  assert(points.size() == projected_points.size());
  assert(image0.channels() == image1.channels());
  //PointSet2 point_set(features);
  KdTree<2, T> kdtree;
  kdtree.Build(Cv2Eigen(features));
  cv::Rect frame0(cv::Point(0, 0), image0.size());
  cv::Rect frame1(cv::Point(0, 0), image1.size());
  cv::Rect r0, r1;
  int num_channels = image0.channels();
  double half_search_range = search_range / 2.0f;

  std::vector<cv::Point_<T>> match_features(points.size(), cv::Point_<T>(-1, -1));
  if (scores != NULL)
    scores->resize(points.size());
  for (size_t i = 0; i < points.size(); ++i)
  {
    //auto in_points =
    //point_set.range_search(to_rect(projected_points[i], search_range));
    std::vector<cv::Point_<T>> in_points;
    std::vector<Eigen::Matrix<T, 2, 1>> eigen_points;
    Eigen::Matrix<T, 2, 1> lower(projected_points[i].x - half_search_range,
                                 projected_points[i].y - half_search_range);
    Eigen::Matrix<T, 2, 1> upper(projected_points[i].x + half_search_range,
                                 projected_points[i].y + half_search_range);
    // Eigen::Matrix<T, 2, 1> lower((int)(projected_points[i].x - half_search_range + 0.5),
    //                              (int)(projected_points[i].y - half_search_range + 0.5));
    // Eigen::Matrix<T, 2, 1> upper = lower + Eigen::Matrix<T, 2, 1>(search_range, search_range);
    kdtree.RangeSearch(lower, upper, std::back_inserter(eigen_points));
    // auto r = to_rect(projected_points[i], search_range);
    // kdtree.RangeSearch(
    //     Eigen::Matrix<T, 2, 1>(r.x, r.y),
    //     Eigen::Matrix<T, 2, 1>(r.x + r.width, r.y + r.height),
    //     std::back_inserter(eigen_points));
    in_points = Eigen2Cv(eigen_points);
    double best_score = psf->worst();
    // if (in_points.begin() == in_points.end())
    // {
    //   std::cout << "Find non match\n";
    // }

    //bool no_found = true;

    for (auto in_it = in_points.begin(); in_it != in_points.end(); ++in_it)
    {
      r0 = to_rect(points[i], patch_size);
      r1 = to_rect(*in_it, patch_size);
      if (r0 <= frame0 && r1 <= frame1)
      {
        double ncc_score = sum((*psf)(image0(r0), image1(r1)));
        if (psf->better(ncc_score, best_score))
        {
          match_features[i] = *in_it;
          best_score = ncc_score;
          //no_found = false;
        }
      }
      // else
      // {
      //   std::cout << r1 << " not in frame" << std::endl;
      // }
    }
    if (scores != NULL)
      (*scores)[i] = best_score / num_channels;
    // if (no_found)
    // {
    //   std::cout << points[i].x << " " << points[i].y << " no found " << in_points.size() << std::endl;
    // }
  }

  return match_features;
}

std::vector<cv::Point2d>
MatchProjectedFeatures(const cv::Mat& image0,
                       const cv::Mat& image1,
                       const std::vector<cv::Point2d> &points,
                       const std::vector<cv::Point2d> &projected_points,
                       const std::vector<cv::Point2d> &features,
                       double corner_strength_weight0,
                       double corner_strength_weight1,
                       int patch_size,
                       int search_range,
                       PatchScoreFunc* psf,
                       std::vector<double>* scores = NULL);

} // furry

#endif // FURRY_COMMON_PATCH
