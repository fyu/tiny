#include "sac_estimator.hpp"
#include <iostream>
#include <memory>
#include <fstream>
#include <opencv2/opencv.hpp>
#include <opencv2/core/eigen.hpp>
#include <ctime>
#include <cstdlib>

using namespace furry;
using namespace Eigen;
using namespace  std;
using namespace cv;

int main(int argc, char** argv)
{
  // //SacEstimator<Line2DSacModel>* sace = new RansacEstimator<Line2DSacModel>();
  // std::shared_ptr<RansacEstimator<Line2DSacModel>> sace(new RansacEstimator<Line2DSacModel>());
  // //auto sace = new RansacEstimator<Line2DSacModel>();
  // sace->residualThreshold(1);
  // vector<Eigen::Vector2d> points;
  // for (int i = 0; i < 1000; ++i)
  // {
  //   Vector2d p(i, i);
  //   p[1] += randn() * 0;
  //   //cout << rand(0, 1) << endl;
  //     //cout << p << endl;
  //   points.push_back(p);
  // }
  // sace->setData(points);
  // for (int i = 0; i < 10; ++i)
  //   sace->run();
  // cout << sace->model() << endl;
  // //delete sace;
  std::srand(std::time(NULL));

  vector<pair<Vector2d, Vector2d>> points;
  ifstream ifs("points_0_0_1_0_nf_10000_r_1_mf_nmf_20000_ws_7_zoom_1.txt");
  string line;
  double ncc;
  while (ifs.good())
  {
    if (ifs.peek() == '#')
      getline(ifs, line);
    else
    {
      pair<Vector2d, Vector2d> pp;
      ifs >> pp.first[0] >> pp.first[1] >> pp.second[0] >> pp.second[1] >> ncc;
      points.push_back(pp);
    }
  }
  points.resize(1000);
  vector<cv::Point2d> points1(points.size());
  vector<cv::Point2d> points2(points.size());
  for (int i = 0; i < points.size(); ++i)
  {
    points1[i].x = points[i].first[0];
    points1[i].y = points[i].first[1];
    points2[i].x = points[i].second[0];
    points2[i].y = points[i].second[1];
  }
  cout << points.size() << " " << points1.size() << endl;
  //Mat cvf = findFundamentalMat(points1, points2, CV_FM_8POINT);
  vector<unsigned char> mask(points.size(), 0);
  Mat cvf = findFundamentalMat(points1, points2, CV_FM_RANSAC, 3, 0.99, mask);
  Matrix3d ef;
  cv::cv2eigen(cvf, ef);
  cout << ef << endl;
  int num_inliers = 0;
  for (int i = 0; i < points.size(); ++i)
  {
    cout << (int)mask[i] << " ";
    num_inliers += (int)mask[i];
  }
  cout << endl;
  cout << num_inliers << endl;

  RansacEstimator<FundmentalMatrixSacModel> ransac;
  ransac.residualThreshold(3);
  ransac.setData(points);
  for (int i = 0; i < 10000; ++i)
  {
    ransac.run();
    if (ransac.numInliers() > 290)
      break;
    //
  }
  cout << "num trials: " << ransac.numRounds() << endl;
  cout << "num inliers: " << ransac.numInliers() << endl;
  cout << "inliers: ";
  auto inlier_mask = ransac.inlierMask();
  for (int i = 0; i < inlier_mask.size(); ++i)
  {
    cout << inlier_mask[i] << " ";
  }
  cout << endl;
  cout << ransac.model() << endl;
  return 0;
}
