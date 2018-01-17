#ifndef FURRY_COMMON_CAMERA_H_
#define FURRY_COMMON_CAMERA_H_

#include <memory>
#include <vector>

#include "opencv2/core/core.hpp"
// #include <gsl/gsl_multiroots.h>
// #include <gsl/gsl_vector.h>
#include <Eigen/Dense>

#include "furry/common/cv.h"
#include "furry/common/pose.h"

namespace furry
{
enum CameraType
	{
		PERSPECTIVE,
		FISH_EYE
	};

struct CameraParam
{
	CameraType type;
	double focal_x;
	double focal_y;
	double center_x;
	double center_y;
	double k1;
	double k2;
	double k3;
	double p1;
	double p2;
	double fov_max;
};

struct GslFuncParam
{
  double k1;
  double k2;
  double k3;
  double p1;
  double p2;

  // distorted position
  double xd;
  double yd;
};

class Camera
{
public:
	Camera(const CameraParam & param);
  virtual ~Camera();

	template <typename T>
	cv::Mat undistort(const cv::Mat & in,
										cv::Size out_size) const
	{
		auto in_size = in.size();
		double x_scale = (double)out_size.width / in_size.width;
		double y_scale = (double)out_size.width / in_size.width;

		//cv::Mat out(out_size, in.type());
		cv::Mat out = cv::Mat::zeros(out_size, in.type());

		for (int x = 0; x < out_size.width; ++x)
		{
			for (int y = 0; y < out_size.height; ++y)
			{

				auto distorted_pos =
					//distortedPos(cv::Point2d(std::atan((x * x_scale - _p.center_x) / _p.focal_x),
					//											 std::atan((y * y_scale - _p.center_y) / _p.focal_y)));
					getDistortedPos(cv::Point2d((x * x_scale - _p.center_x) / _p.focal_x,
																	 (y * y_scale - _p.center_y) / _p.focal_y));
        distorted_pos.x = _p.focal_x * distorted_pos.x + _p.center_x;
        distorted_pos.y = _p.focal_y * distorted_pos.y + _p.center_y;
				out.at<T>(cv::Point(x, y)) = m_at_bilinear<T>(in, distorted_pos);
			}
		}
		return out;
	}
  //cv::Point2d getUndistortedPos(cv::Point2d);
  virtual cv::Point2d getDistortedPos(cv::Point2d) const = 0;
  virtual Eigen::Vector2d getDistortedPos(Eigen::Vector2d) const = 0;
  virtual void print_info() const = 0;

  Eigen::Matrix3d intrinsics() const;

  // Project the point to z = 1 plane with center (0, 0)
  cv::Point2d normalize(const cv::Point2d&) const;
  cv::Point2d denormalize(const cv::Point2d&) const;
  Eigen::Vector2d normalize(const Eigen::Vector2d&) const;
  Eigen::Vector2d denormalize(const Eigen::Vector2d&) const;

  Eigen::Vector3d PixelInWorld(const Pose &camera_pose,
                               const Eigen::Vector2d &p,
                               bool distorted = false);
  double focal() const;

protected:
	CameraParam _p;

  // const gsl_multiroot_fdfsolver_type* _gsl_solver_type;
  // gsl_multiroot_fdfsolver* _gsl_solver;
  // gsl_multiroot_function_fdf _gsl_func;
  // gsl_vector* _gsl_init;
  // GslFuncParam _gsl_param;
};

class PerspectiveCamera: public Camera
{
public:
	PerspectiveCamera(const CameraParam & param);
  virtual void print_info() const;
	virtual cv::Point2d getDistortedPos(cv::Point2d) const;
  virtual Eigen::Vector2d getDistortedPos(Eigen::Vector2d) const;
};

class FishEyeCamera : public Camera
{
public:
	FishEyeCamera(const CameraParam & param) : Camera(param) {}
  virtual void print_info() const;
	virtual cv::Point2d getDistortedPos(cv::Point2d) const { return cv::Point2d(0, 0); }
  virtual Eigen::Vector2d getDistortedPos(Eigen::Vector2d) const { return Eigen::Vector2d(0, 0); }
};

template <template <typename> class Ptr>
Ptr<Camera> make_camera(const CameraParam & param)
{
  switch (param.type)
	{
	case PERSPECTIVE:
		return Ptr<Camera>(new PerspectiveCamera(param));
	case FISH_EYE:
		return Ptr<Camera>(new FishEyeCamera(param));
	}
  // To mute the compiler warning
  return Ptr<Camera>(new PerspectiveCamera(param));
}

////////////////////////////////////////////////////////////
// Camera Inline
////////////////////////////////////////////////////////////

inline cv::Point2d
Camera::normalize(const cv::Point2d& p) const
{
  return cv::Point2d((p.x - _p.center_x) / _p.focal_x,
                     (p.y - _p.center_y) / _p.focal_y);
}

inline cv::Point2d
Camera::denormalize(const cv::Point2d& p) const
{
  return cv::Point2d(p.x * _p.focal_x + _p.center_x,
                     p.y * _p.focal_y + _p.center_y);
}

inline Eigen::Vector2d
Camera::normalize(const Eigen::Vector2d& p) const
{
  return Eigen::Vector2d((p.x() - _p.center_x) / _p.focal_x,
                         (p.y() - _p.center_y) / _p.focal_y);
}

inline Eigen::Vector2d
Camera::denormalize(const Eigen::Vector2d& p) const
{
  return Eigen::Vector2d(p.x() * _p.focal_x + _p.center_x,
                         p.y() * _p.focal_y + _p.center_y);
}

inline Eigen::Matrix3d
Camera::intrinsics() const
{
  Eigen::Matrix3d m = Eigen::Matrix3d::Identity();
  m(0, 0) = _p.focal_x;
  m(1, 1) = _p.focal_y;
  m(0, 2) = _p.center_x;
  m(1, 2) = _p.center_y;

  return m;
}

} // furry

#endif // FURRY_COMMON_CAMERA
