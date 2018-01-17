#include "furry/common/geometry.h"

//#include <iostream>

namespace furry {

bool LineEndPoints(const Eigen::Vector3d &line, int width, int height,
                   Eigen::Vector2d *p0, Eigen::Vector2d *p1) {
  Eigen::Vector2d c0 =
      line.cross(Eigen::Vector3d(-height, 0, 0)).hnormalized();
  Eigen::Vector2d c1 =
      line.cross(Eigen::Vector3d(0, width, 0)).hnormalized();
  Eigen::Vector2d c2 =
      line.cross(Eigen::Vector3d(-height, 0, height * width)).hnormalized();
  Eigen::Vector2d c3 =
      line.cross(Eigen::Vector3d(0, width, -height * width)).hnormalized();

  //std::cout << c0 << '\n' <<  c1 << '\n' << c2 << '\n' <<  c3 << std::endl;

  bool found = false;
  if (c0.y() >= 0 && c0.y() <= height) {
    found = true;
    *p0 = c0;
  }
  if (c1.x() >= 0 && c1.x() <= width) {
    if (!found) *p0 = c1;
    else *p1 = c1;
    found = true;
  }
  if (c2.y() >= 0 && c2.y() <= height) {
    if (!found) *p0 = c2;
    else *p1 = c2;
    found = true;
  }
  if (c3.x() >= 0 && c3.x() <= width) {
    if (!found) *p0 = c3;
    else *p1 = c3;
    found = true;
  }
  return found;
}

Eigen::Matrix3d Quaternion2RotationMatrix(double x, double y, double z, double w)
{
  Eigen::Matrix3d m = Eigen::Matrix3d::Identity();
  m(0, 0) = 1 - 2 * (y * y + z * z);
  m(0, 1) = 2 * (x * y - z * w);
  m(0, 2) = 2 * (x * z + y * w);
  m(1, 0) = 2 * (x * y + z * w);
  m(1, 1) = 1 - 2 * (x * x + z * z);
  m(1, 2) = 2 * (y * z - x * w);
  m(2, 0) = 2 * (x * z - y * w);
  m(2, 1) = 2 * (y * z + x * w);
  m(2, 2) = 1 - 2 * (x * x + y * y);

  return m;
}

} // furry
