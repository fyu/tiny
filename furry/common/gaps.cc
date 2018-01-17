#include "furry/common/gaps.h"

namespace furry {

cv::Mat R2Image2CvImage(const R2Image *image) {
  auto flipped_image =  cv::Mat(image->Height(),
                                image->Width(),
                                CV_MAKE_TYPE(CV_8U, image->NComponents()),
                                (void*)image->Pixels(),
                                image->RowSize());
  cv::Mat cvimage;
  cv::flip(flipped_image, cvimage, 0);
  cv::cvtColor(cvimage, cvimage, CV_RGB2BGR);
  return cvimage;
}

////////////////////////////////////////////////////////////
// R3Shapes
////////////////////////////////////////////////////////////

R3Quaternion
Conjugate(const R3Quaternion &q) {
  return R3Quaternion(q.A(), -q.B(), -q.C(), -q.D());
}

R3Quaternion
RotateBetween(const R3Quaternion &from,
              const R3Quaternion &to) {
  return to * Conjugate(from);
}

R3Vector
Rotate(const R3Quaternion &q, const R3Vector &v) {
  auto m = q.Matrix();
  return m * v;
}

R3Point
Rotate(const R3Quaternion &q, const R3Point &p) {
  auto m = q.Matrix();
  return m * p;
}

Eigen::Vector3d R3QuaternionToAngleAxis(const R3Quaternion &q) {
  Eigen::AngleAxisd aa(To<Eigen::Quaterniond>(q));
  return aa.angle() * aa.axis();
}

double CalcDistance2(const R3Point &p0, const R3Point &p1) {
  return
      (p0.X() - p1.X()) * (p0.X() - p1.X()) +
      (p0.Y() - p1.Y()) * (p0.Y() - p1.Y()) +
      (p0.Z() - p1.Z()) * (p0.Z() - p1.Z());
}

double CalcDistance(const R3Point &p0, const R3Point &p1) {
  return sqrt(CalcDistance2(p0, p1));
}

} // furry

////////////////////////////////////////////////////////////
// IO
////////////////////////////////////////////////////////////

std::ostream& operator << (std::ostream &os, const R3Plane &p) {
  os << p.A() << ' ' << p.B() << ' ' << p.C() << ' ' << p.D();
  return os;
}

std::ostream& operator << (std::ostream &os, const R3Point &p) {
  os << p.X() << ' ' << p.Y() << ' ' << p.Z();
  return os;
}

std::ostream& operator << (std::ostream &os, const R2Point &p) {
  os << p.X() << ' ' << p.Y();
  return os;
}

std::ostream& operator << (std::ostream &os, const R3Quaternion &q) {
  os << q.A() << ' ' << q.B() << ' ' << q.C() << ' ' << q.D();
  return os;
}

std::ostream& operator << (std::ostream &os, const R3Vector &v) {
  os << v.X() << ' ' << v.Y() << ' ' << v.Z();
  return os;
}
