#include <cmath>
#include <limits>
#include <fstream>
#include <string>
#include <cstdio>

#include "furry/common/camera.h"
#include "furry/common/cv.h"

using std::fstream;
//using std::shared_ptr;
using std::string;
using std::vector;
using std::printf;
using cv::Point2d;

namespace furry
{

Camera::Camera(const CameraParam & param): _p(param)
{
  // _gsl_solver_type = gsl_multiroot_fdfsolver_hybridsj;
  // _gsl_solver = gsl_multiroot_fdfsolver_alloc(_gsl_solver_type, 2);
  // _gsl_init = gsl_vector_alloc(2);

  // _gsl_param.k1 = _p.k1;
  // _gsl_param.k2 = _p.k2;
  // _gsl_param.k3 = _p.k3;
  // _gsl_param.p1 = _p.p1;
  // _gsl_param.p2 = _p.p2;
  // _gsl_func.n = 2;
  // _gsl_func.params = &_gsl_param;
}

Camera::~Camera()
{
  // gsl_multiroot_fdfsolver_free(_gsl_solver);
  // gsl_vector_free(_gsl_init);
}

double Camera::focal() const {
  return (_p.focal_x + _p.focal_y) / 2;
}

Eigen::Vector3d
Camera::PixelInWorld(const Pose &camera_pose,
                     const Eigen::Vector2d &p,
                     bool distorted) {
  if (distorted) {
    std::cerr << "Error: Undistorting pixel is not supported\n";
  }

  auto towards = camera_pose.towards();
  auto up = camera_pose.up();
  auto right = towards.cross(up);

  auto np = normalize(p);
  // double dx = 2 * (p.x() - 968) / 1936;
  // double dy = 2 * (p.y() - 1296) / 2592;
  Eigen::Vector3d world_position = camera_pose.position();
  world_position += np.x() * right;
  world_position += np.y() * up;

  return world_position;
}
 
// Point2d
// Camera::getUndistortedPos(Point2d point)
// {
//   _gsl_param.xd = point.x;
//   _gsl_param.yd = point.y;
//   gsl_vector_set(_gsl_init, 0, point.x);
//   gsl_vector_set(_gsl_init, 1, point.y);
//   gsl_multiroot_fdfsolver_set(_gsl_solver, &_gsl_func, _gsl_init);

//   int status, iter = 0;
//   do
//   {
//     iter++;
//     status = gsl_multiroot_fdfsolver_iterate(_gsl_solver);
//     if (status)
//       break;
//     status = gsl_multiroot_test_residual(_gsl_solver->f, 1e-8);
//   }
//   while(status == GSL_CONTINUE && iter < 1000);

//   Point2d result;
//   result.x = gsl_vector_get(_gsl_solver->x, 0);
//   result.y = gsl_vector_get(_gsl_solver->x, 1);
//   return result;
// }

// int
// perspective_distort_f (const gsl_vector* gsl_xy, void* params, gsl_vector* f)
// {
//   GslFuncParam* p = (GslFuncParam*)params;
//   const double x = gsl_vector_get(gsl_xy, 0);
//   const double y = gsl_vector_get(gsl_xy, 1);

//   double x2 = x * x;
//   double y2 = y * y;
//   double r2 = x2 + y2;
//   double xy = x * y;

//   double radial_scale = 1 + p->k1 * r2 + p->k2 * r2 * r2 + p->k3 * r2 * r2 * r2;
//   gsl_vector_set(f, 0,
//                  radial_scale * x + 2 * p->p1 * xy + p->p2 * (r2 + 2 * x2) - p->xd);
//   gsl_vector_set(f, 1,
//                  radial_scale * y + p->p1 * (r2 + 2 * y2) + 2 * p->p2 * xy - p->yd);
//   return GSL_SUCCESS;
// }

// int
// perspective_distort_df(const gsl_vector* gsl_xy, void* params, gsl_matrix* J)
// {
//   GslFuncParam* p = (GslFuncParam*)params;
//   const double x = gsl_vector_get(gsl_xy, 0);
//   const double y = gsl_vector_get(gsl_xy, 1);

//   double x2 = x * x;
//   double y2 = y * y;
//   double r2 = x2 + y2;
//   double xy = x * y;

//   double radial_f = 1 + p->k1 * r2 + p->k2 * r2 * r2 + p->k3 * r2 * r2 * r2;
//   double radial_dfx = p->k1 * 2 * x + p->k2 * 4 * r2 * x + p->k3 * 6 * r2 * r2 * x;
//   double radial_dfy = p->k1 * 2 * y + p->k2 * 4 * r2 * y + p->k3 * 6 * r2 * r2 * y;
//   gsl_matrix_set(J, 0, 0,
//                  radial_f + x * radial_dfx + 2 * p->p1 * y + p->p2 * 6 * x);
//   gsl_matrix_set(J, 0, 1,
//                  radial_dfy * x + 2 * p->p1 * x + p->p2 * 2 * y);
//   gsl_matrix_set(J, 1, 0,
//                  radial_dfx * y + p->p1 * 2 * x + 2 * p->p2 * y);
//   gsl_matrix_set(J, 1, 1,
//                  radial_f + radial_dfy * y + p->p1 * 6 * y + 2 * p->p2 * x);
//   return GSL_SUCCESS;
// }

// int
// perspective_distort_fdf(const gsl_vector* gsl_xy,
//                                void* params,
//                                gsl_vector* f,
//                                gsl_matrix* J)
// {
//   GslFuncParam* p = (GslFuncParam*)params;
//   const double x = gsl_vector_get(gsl_xy, 0);
//   const double y = gsl_vector_get(gsl_xy, 1);

//   double x2 = x * x;
//   double y2 = y * y;
//   double r2 = x2 + y2;
//   double xy = x * y;

//   double radial_f = 1 + p->k1 * r2 + p->k2 * r2 * r2 + p->k3 * r2 * r2 * r2;
//   double radial_dfx = p->k1 * 2 * x + p->k2 * 4 * r2 * x + p->k3 * 6 * r2 * r2 * x;
//   double radial_dfy = p->k1 * 2 * y + p->k2 * 4 * r2 * y + p->k3 * 6 * r2 * r2 * y;

//   gsl_vector_set(f, 0,
//                  radial_f * x + 2 * p->p1 * xy + p->p2 * (r2 + 2 * x2) - p->xd);
//   gsl_vector_set(f, 1,
//                  radial_f * y + p->p1 * (r2 + 2 * y2) + 2 * p->p2 * xy - p->yd);

//   gsl_matrix_set(J, 0, 0,
//                  radial_f + x * radial_dfx + 2 * p->p1 * y + p->p2 * 6 * x);
//   gsl_matrix_set(J, 0, 1,
//                  radial_dfy * x + 2 * p->p1 * x + p->p2 * 2 * y);
//   gsl_matrix_set(J, 1, 0,
//                  radial_dfx * y + p->p1 * 2 * x + 2 * p->p2 * y);
//   gsl_matrix_set(J, 1, 1,
//                  radial_f + radial_dfy * y + p->p1 * 6 * y + 2 * p->p2 * x);
//   return GSL_SUCCESS;
// }

PerspectiveCamera::PerspectiveCamera(const CameraParam& param): Camera(param)
{
  // _gsl_func.f = perspective_distort_f;
  // _gsl_func.df = perspective_distort_df;
  // _gsl_func.fdf = perspective_distort_fdf;
}

cv::Point2d
PerspectiveCamera::getDistortedPos(cv::Point2d p) const
{
	//double r_max = tan( _p.fov_max / 360.0 * M_PI);
	//double r2_max = r_max * r_max;

	double x2 = p.x * p.x;
	double y2 = p.y * p.y;
	double r2 = x2 + y2;

	//if (r2 > r2_max) return nan_point;

	double xy = p.x * p.y;

	double radial_scale = 1 + _p.k1 * r2 + _p.k2 * r2 * r2 + _p.k3 * r2 * r2 * r2;
	double x_distorted = radial_scale * p.x + 2 * _p.p1 * xy + _p.p2 * (r2 + 2 * x2);
	double y_distorted = radial_scale * p.y + _p.p1 * (r2 + 2 * y2) + 2 * _p.p2 * xy;

	//x_distorted = _p.focal_x * x_distorted + _p.center_x;
	//y_distorted = _p.focal_y * y_distorted + _p.center_y;

	return cv::Point2f(x_distorted, y_distorted);
}

Eigen::Vector2d
PerspectiveCamera::getDistortedPos(Eigen::Vector2d p) const
{
  //double r_max = tan( _p.fov_max / 360.0 * M_PI);
	//double r2_max = r_max * r_max;
  double px = p.x();
  double py = p.y();

	double x2 = px * px;
	double y2 = py * py;
	double r2 = x2 + y2;

	//if (r2 > r2_max) return nan_point;

	double xy = px * py;

	double radial_scale = 1 + _p.k1 * r2 + _p.k2 * r2 * r2 + _p.k3 * r2 * r2 * r2;
	double x_distorted = radial_scale * px + 2 * _p.p1 * xy + _p.p2 * (r2 + 2 * x2);
	double y_distorted = radial_scale * py + _p.p1 * (r2 + 2 * y2) + 2 * _p.p2 * xy;

  return Eigen::Vector2d(x_distorted, y_distorted);
}

void PerspectiveCamera::print_info() const
{
  printf("=== Camera Info ===\n\n");
  printf("distortion type: Perspective\n");
  printf("focal x: %f y: %f\n", _p.focal_x, _p.focal_y);
  printf("center x: %f y: %f\n", _p.center_x, _p.center_y);
  printf("k1: %f k2: %f k3: %f\n", _p.k1, _p.k2, _p.k3);
  printf("p1: %f p2: %f\n", _p.p1, _p.p2);
  printf("max fov: %f\n\n", _p.fov_max);
}

void FishEyeCamera::print_info() const
{
}

}

