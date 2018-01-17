#include "furry/common/patch.h"
#include "furry/common/cv.h"
#include "furry/common/point_set.h"

#include <iostream>
#include <limits>
#include <algorithm>
#include <cmath>

namespace furry
{

cv::Rect
to_rect(double x, double y, double width, double height)
{
  return cv::Rect(x - width / 2.0 + 0.5,
                  y - height / 2.0 + 0.5,
                  width + 0.5,
                  height + 0.5);
}

cv::Rect
to_rect(const cv::Point2d& p, double s)
{
  return to_rect(p.x, p.y, s, s);
}


cv::Mat
make_patch(const cv::Mat& mat, double x, double y, double width, double height)
{
  return mat(to_rect(x, y, width, height));
}

cv::Scalar
ssd(const cv::Mat& p0, const cv::Mat& p1)
{
  cv::Scalar s = cv::sum((p0 - p1).mul(p0 - p1));
  int area = p0.size().width * p0.size().height;
  return cv::Scalar(s[0] / area, s[1] / area, s[2] / area, s[3] / area);
}

cv::Scalar
cc(const cv::Mat& p0, const cv::Mat& p1)
{
  return cv::sum(p0.mul(p1));
}

cv::Scalar
ncc(const cv::Mat& p0, const cv::Mat& p1)
{
  cv::Mat d0 = p0 - cv::mean(p0);
  cv::Mat d1 = p1 - cv::mean(p1);
  cv::Scalar cov = cv::sum(d0.mul(d1));
  cv::Scalar sigma;
  cv::sqrt(cc(d0, d0).mul(cc(d1, d1)), sigma);
  return cov / sigma;
}

cv::Scalar
nccm(const cv::Mat& p0, const cv::Mat& p1)
{
  cv::Scalar mu0 = cv::mean(p0);
  cv::Scalar mu1 = cv::mean(p1);
  cv::Mat d0 = p0 - mu0;
  cv::Mat d1 = p1 - mu1;
  cv::Scalar cov = cv::sum(d0.mul(d1));
  cv::Scalar sigma;
  cv::sqrt(cc(d0, d0).mul(cc(d1, d1)), sigma);
  return cov / sigma * (cv::Scalar(1, 1, 1, 1) -
                        furry::abs(mu0 - mu1) / (mu0 + mu1));
}

cv::Scalar
PatchScoreFunc::operator () (const cv::Mat& p0, const cv::Mat& p1) const
{
  cv::Mat p0p = p0;
  cv::Mat p1p = p1;
  if (p0.depth() != CV_64F)
    p0.convertTo(p0p, CV_64F);
  if (p1.depth() != CV_64F)
    p1.convertTo(p1p, CV_64F);
  return compute(p0p, p1p);
}

PatchScoreFunc::~PatchScoreFunc()
{
}

cv::Scalar
PatchSsd::compute(const cv::Mat& p0, const cv::Mat& p1) const
{
  return ssd(p0, p1);
}

bool
PatchSsd::better(double s0, double s1) const
{
  return s0 < s1;
}

double
PatchSsd::worst() const
{
  return std::numeric_limits<double>::max();
}

cv::Scalar
PatchNcc::compute(const cv::Mat& p0, const cv::Mat& p1) const
{
  return ncc(p0, p1);
}

bool
PatchNcc::better(double s0, double s1) const
{
  return s0 > s1;
}

double
PatchNcc::worst() const
{
  //return std::numeric_limits<double>::min();
  return -1e5;
}


cv::Scalar
PatchNccm::compute(const cv::Mat& p0, const cv::Mat& p1) const
{
  return nccm(p0, p1);
}

bool
PatchNccm::better(double s0, double s1) const
{
  return s0 > s1;
}

double
PatchNccm::worst() const
{
  //return std::numeric_limits<double>::min();
  return -1e5;
}

cv::Point2d
adjust_patch_match(const cv::Mat& image0,
                      const cv::Mat& image1,
                      const cv::Point2d& p0,
                      const cv::Point2d& p1,
                      int patch_size,
                      int range,
                      PatchScoreFunc* psf,
                      double* final_score)
{
  //assert(image0.depth() == image1.depth() == CV_64F);
  cv::Mat image0d = image0;
  cv::Mat image1d = image1;
  if (image0.depth() != CV_64F)
    image0.convertTo(image0d, CV_64F);
  if (image1.depth() != CV_64F)
    image1.convertTo(image1d, CV_64F);

  cv::Rect r0 = to_rect(p0, patch_size);
  cv::Rect r1 = to_rect(p1, patch_size);
  cv::Rect frame0 = imframe(image0);
  cv::Rect frame1 = imframe(image1);

  if (!(r0 <= frame0))
    return p1;

  cv::Mat patch0 = image0d(r0);

  cv::Rect best_rect = r1;
  //Scalar score = score_func(patch0, image1(r1));
  cv::Scalar score = (*psf)(patch0, image1d(r1));
  double best_score = score[0] + score[1] + score[2];

  cv::Rect crect;
  for (int dx = -range; dx <= range; ++dx)
  {
    crect = r1 + cv::Point(dx, 0);
    if (crect <= frame1)
    {
      for (int dy = -range; dy <= range; ++dy)
      {
        crect = r1 + cv::Point(dx, dy);
        if (crect <= frame1)
        {
          score = (*psf)(patch0, image1d(crect));
          double total = score[0] + score[1] + score[2];
          if (psf->better(total, best_score))
          {
            best_rect = crect;
            best_score = total;
          }
        }
      }
    }
  }

  if (final_score != NULL) *final_score = best_score / 3;

  return cv::Point2d(best_rect.x + patch_size / 2.0,
                     best_rect.y + patch_size / 2.0);
}

// std::vector<cv::Point2d>
// MatchProjectedFeatures(const cv::Mat& image0,
//                        const cv::Mat& image1,
//                        const std::vector<cv::Point2d> &points,
//                        const std::vector<cv::Point2d> &projected_points,
//                        const std::vector<cv::Point2d> &features,
//                        int patch_size,
//                        int search_range,
//                        PatchScoreFunc* psf,
//                        std::vector<double>* scores)
// {
//   assert(points.size() == projected_points.size());
//   assert(image0.channels() == image1.channels());
//   PointSet2 point_set(features);
//   cv::Rect frame0(cv::Point(0, 0), image0.size());
//   cv::Rect frame1(cv::Point(0, 0), image1.size());
//   cv::Rect r0, r1;
//   int num_channels = image0.channels();

//   std::vector<cv::Point2d> match_features(points.size(), cv::Point2d(-1, -1));
//   if (scores != NULL)
//     scores->resize(points.size());
//   for (size_t i = 0; i < points.size(); ++i)
//   {
//     auto in_points =
//       point_set.range_search(to_rect(projected_points[i], search_range));
//     double best_score = psf->worst();
//     // if (in_points.begin() == in_points.end())
//     // {
//     //   std::cout << "Find non match\n";
//     // }

//     //bool no_found = true;

//     for (auto in_it = in_points.begin(); in_it != in_points.end(); ++in_it)
//     {
//       r0 = to_rect(points[i], patch_size);
//       r1 = to_rect(*in_it, patch_size);
//       if (r0 <= frame0 && r1 <= frame1)
//       {
//         double ncc_score = sum((*psf)(image0(r0), image1(r1)));
//         if (psf->better(ncc_score, best_score))
//         {
//           match_features[i] = *in_it;
//           best_score = ncc_score;
//           //no_found = false;
//         }
//       }
//       // else
//       // {
//       //   std::cout << r1 << " not in frame" << std::endl;
//       // }
//     }
//     if (scores != NULL)
//       (*scores)[i] = best_score / num_channels;
//     // if (no_found)
//     // {
//     //   std::cout << points[i].x << " " << points[i].y << " no found " << in_points.size() << std::endl;
//     // }
//   }

//   return match_features;
// }

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
                       std::vector<double>* scores)
{
  //printf("num of features: %lu\n", features.size());
  assert(points.size() == projected_points.size());
  assert(image0.channels() == image1.channels());
  PointSet2 point_set(features);
  cv::Rect frame0(cv::Point(0, 0), image0.size());
  cv::Rect frame1(cv::Point(0, 0), image1.size());
  cv::Rect r0, r1;
  int num_channels = image0.channels();

  cv::Mat eigen0, eigen1;
  cv::cornerEigenValsAndVecs(image0, eigen0, 3, 3);
  cv::cornerEigenValsAndVecs(image1, eigen1, 3, 3);

  std::vector<cv::Point2d> match_features(points.size(), cv::Point2d(-1, -1));
  if (scores != NULL)
    scores->resize(points.size());
  for (size_t i = 0; i < points.size(); ++i)
  {
    auto in_points =
      point_set.range_search(to_rect(projected_points[i], search_range));
    //double best_score = psf->worst();
    double best_score = std::numeric_limits<double>::max();
    // if (in_points.begin() == in_points.end())
    // {
    //   std::cout << "Find non match\n";
    // }

    //bool no_found = true;

    //auto e00 = eigen0.at<cv::Vec6f>(points[i])[0];
    auto e01 = eigen0.at<cv::Vec6f>(points[i])[1];

    // printf("eigen values %f %f %f %f\n", eigen_value00, eigen_value01,
    //        corner_strength_weight0,
    //        corner_strength_weight1);

    for (auto in_it = in_points.begin(); in_it != in_points.end(); ++in_it)
    {
      r0 = to_rect(points[i], patch_size);
      r1 = to_rect(*in_it, patch_size);
      // auto e10 = eigen1.at<cv::Vec6f>(*in_it)[0];
      auto e11 = eigen1.at<cv::Vec6f>(*in_it)[1];
      if (r0 <= frame0 && r1 <= frame1)
      {
        // double strength_score =
        //     corner_strength_weight0 *
        //     (eigen_value00 - eigen_value10) * (eigen_value00 - eigen_value10) +
        //     corner_strength_weight1 *
        //     (eigen_value01 - eigen_value11) * (eigen_value01 - eigen_value11);
        //printf("strength score: %f\n", strength_score);
        double strength_score = pow((e01 - e11) * (e01 - e11),
                                    corner_strength_weight0);
#ifndef MAX
#define MAX(a, b) ((a) > (b))? (a) : (b);
#endif
        double ncc_score = //1 - sum((*psf)(image0(r0), image1(r1))) / 3;
            std::exp(1 - MAX(sum((*psf)(image0(r0), image1(r1))) / 3, 0)) * strength_score;
        //printf("total score: %f\n", ncc_score);
        //if (!psf->better(ncc_score, best_score))
        if (ncc_score < best_score)
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

} // furry
