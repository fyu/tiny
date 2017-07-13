#include "furry/common/cv.h"

#include <iterator>
#include <limits>

#include <glog/logging.h>

#include "furry/common/memory.h"
#include "furry/common/point.h"

namespace furry
{

const cv::Scalar kCvWhite = CV_RGB(255, 255, 255);
const cv::Scalar kCvBlack = CV_RGB(0, 0, 0);
const cv::Scalar kCvBlue = CV_RGB(0, 0, 255);
const cv::Scalar kCvRed = CV_RGB(255, 0, 0);
const cv::Scalar kCvGreen = CV_RGB(0, 255, 0);
const cv::Scalar kCvMagenta = CV_RGB(255, 0, 255);
const cv::Scalar kCvCyan = CV_RGB(0, 255, 255);
const cv::Scalar kCvYellow = CV_RGB(255, 255, 0);
const cv::Scalar kCvGray = CV_RGB(128, 128, 128);
const cv::Scalar kCvLightGray = CV_RGB(211, 211, 211);

extern const int kCvColorPlatteSize = 6;
extern const cv::Scalar kCvColorPlatte[kCvColorPlatteSize] = {
  kCvRed, kCvGreen, kCvBlue, kCvYellow, kCvMagenta, kCvCyan};

////////////////////////////////////////////////////////////
// Vec
////////////////////////////////////////////////////////////

double RgbToLuminosity(const cv::Vec3b &c) {
  return c[0] * 0.21 + c[1] * 0.72 + c[2] * 0.07;
}

double BgrToLuminosity(const cv::Vec3b &c) {
  return c[2] * 0.21 + c[1] * 0.72 + c[0] * 0.07;
}

////////////////////////////////////////////////////////////
// Point_
////////////////////////////////////////////////////////////

template <> bool
equal<cv::Point2d>(const cv::Point2d& p0, const cv::Point2d& p1)
{
  return equal(p0.x, p1.x) && equal(p0.y, p1.y);
}

std::ostream& operator << (std::ostream& os, const cv::Point2d& p)
{
  os << p.x << ' ' << p.y;
  return os;
}

std::ostream& operator << (std::ostream& os, const cv::Point3d& p)
{
  os << p.x << ' ' << p.y << ' ' << p.z;
  return os;
}

cv::KeyPoint FlipY(const cv::KeyPoint &p, double height) {
  cv::KeyPoint pp(p);
  pp.pt.y = height - pp.pt.y;
  return pp;
}

////////////////////////////////////////////////////////////
// Feature 2D
////////////////////////////////////////////////////////////

std::ostream& operator << (std::ostream &os, const cv::KeyPoint &p) {
  os << p.pt.x << ' ' << p.pt.y << ' ' << p.size << ' ' << p.angle << ' '
     << p.response << ' ' << p.octave << ' ' << p.class_id;
  return os;
}

////////////////////////////////////////////////////////////
// Point_
////////////////////////////////////////////////////////////

#if 0

// A reimplementation of WritePoints

int WritePoints(const std::string &filename,
                const std::vector<cv::Point2f> &points) {
  static const size_t kBufferSize = 1024;
  static const size_t kElemSize = sizeof(float);
  static const size_t kPointSize = 2 * kElemSize;
  static const size_t kNumBytes = kBufferSize * kPointSize;

  auto buffer = ScopedArray<float>(new float[kBufferSize * 2]);

  int file = open(filename.c_str(),
                  O_WRONLY | O_CREAT | O_TRUNC,
                  S_IRUSR | S_IWUSR);
  if (file < 0) {
    return 1;
  }

  auto first = points.begin();
  auto last = points.end();
  while (first != last) {
    size_t float_to_write = 0;
    unsigned char *char_buffer = reinterpret_cast<unsigned char*>(buffer.Get());
    while ((float_to_write < kBufferSize * 2) && (first != last)) {
      buffer[float_to_write] = first->x;
      buffer[float_to_write + 1] = first->y;
      float_to_write += 2;
    }
    size_t bytes_to_write = float_to_write * kElemSize;
    size_t bytes_written;
    do {
      bytes_written = write(file, char_buffer, bytes_to_write);
      bytes_to_write -= bytes_written;
      char_buffer += bytes_written;
    } while (bytes_to_write);
  }

  close(file);

  return 0;
}

#endif

////////////////////////////////////////////////////////////
// Scalar
////////////////////////////////////////////////////////////

cv::Scalar
operator / (const cv::Scalar& s0, const cv::Scalar& s1)
{
  cv::Scalar result(0, 0, 0, 0);
  for (int i = 0; i < 4; ++i)
  {
    if (fabs(s1[i]) > std::numeric_limits<double>::epsilon())
    {
      result[i] = s0[i] / s1[i];
    }
  }
  return result;
}

cv::Scalar
operator / (const cv::Scalar &s, double d) {
  return cv::Scalar(s[0] / d, s[1] / d, s[2] / d, s[3] / d);
}

cv::Scalar
operator * (const cv::Scalar& s0, const cv::Scalar& s1)
{
  return cv::Scalar(s0[0] * s1[0], s0[1] * s1[1], s0[2] * s1[2], s0[3] * s1[3]);
}

cv::Scalar
operator + (const cv::Scalar& s0, const cv::Scalar& s1)
{
  return cv::Scalar(s0[0] + s1[0], s0[1] + s1[1], s0[2] + s1[2], s0[3] + s1[3]);
}

cv::Scalar
operator - (const cv::Scalar& s0, const cv::Scalar& s1)
{
  return cv::Scalar(s0[0] - s1[0], s0[1] - s1[1], s0[2] - s1[2], s0[3] - s1[3]);
}

cv::Scalar
abs(const cv::Scalar& s)
{
  return cv::Scalar(fabs(s[0]), fabs(s[1]), fabs(s[2]), fabs(s[3]));
}

double
sum(const cv::Scalar& s)
{
  return s[0] + s[1] + s[2] + s[3];
}

////////////////////////////////////////////////////////////
// Size
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// Mat
////////////////////////////////////////////////////////////

bool WriteImage(const cv::Mat &image, FILE *file) {
  int type = image.type();
  int rows = image.rows;
  int cols = image.cols;
  // LOG(INFO) << "type: " << type << " rows: " << rows << " cols: " << cols;
  fwrite(&type, sizeof(type), 1, file);
  fwrite(&rows, sizeof(rows), 1, file);
  fwrite(&cols, sizeof(cols), 1, file);
  CHECK(image.isContinuous());
  fwrite(image.ptr(), image.elemSize(), image.total(), file);
  return true;
}

bool ReadImage(FILE *file, cv::Mat *image) {
  int type, rows, cols;
  fread(&type, sizeof(type), 1, file);
  fread(&rows, sizeof(rows), 1, file);
  fread(&cols, sizeof(cols), 1, file);
  // LOG(INFO) << "type: " << type << " rows: " << rows << " cols: " << cols;
  image->create(rows, cols, type);
  CHECK(image->isContinuous());
  fread(image->ptr(), image->elemSize(), image->total(), file);
  return true;
}

cv::Rect
imframe(const cv::Mat& image)
{
  return cv::Rect(0, 0, image.size().width, image.size().height);
}

////////////////////////////////////////////////////////////
// Drawing functions
////////////////////////////////////////////////////////////

cv::Mat
DrawImages(cv::InputArray image0,
           cv::InputArray image1,
           cv::Size image_size) {
  assert(image0.type() == image1.type());
  auto type = image0.type();
  cv::Size size;
  cv::Mat copy0, copy1;
  if (image_size.width < 0 || image_size.height < 0) {
    assert(image0.size() == image1.size());
    size = image0.size();
    copy0 = image0.getMat();
    copy1 = image1.getMat();
  } else {
    size = image_size;
    cv::resize(image0, copy0, size);
    cv::resize(image1, copy1, size);
  }
  cv::Mat pair(size.height, size.width * 2, type);
  copy0.copyTo(pair(cv::Rect(0, 0, size.width, size.height)));
  copy1.copyTo(pair(cv::Rect(size.width, 0, size.width, size.height)));
  return pair;
}

void DrawMatches(const std::vector<cv::Point2f> &points0,
                 const std::vector<cv::Point2f> &points1,
                 const std::vector<cv::Scalar> &colors0,
                 const std::vector<cv::Scalar> &colors1,
                 int point_radius,
                 int line_width,
                 cv::Mat *image) {
  std::function<cv::Scalar (int)> pick_color0, pick_color1;;
  if (colors0.size() == 1)
    pick_color0 = [&colors0](int i) { return colors0[0]; };
  else
    pick_color0 = [&colors0](int i) { return colors0[i]; };
  if (colors1.size() == 0)
    pick_color1 = pick_color0;
  else if (colors1.size() == 1)
    pick_color1 = [&colors1](int i) { return colors1[0]; };
  else
    pick_color1 = [&colors1](int i) { return colors1[i]; };
  for (size_t i = 0; i < points0.size(); ++i)
    circle(*image, points0[i], point_radius, pick_color0(i), -1);
  for (size_t i = 0; i < points1.size(); ++i)
    circle(*image, points1[i], point_radius, pick_color1(i), -1);
  if (line_width > 0)
    for (size_t i = 0; i < points0.size(); ++i)
      line(*image, points0[i], points1[i], pick_color0(i), line_width);
}

cv::Scalar RandomColor() {
  return cv::Scalar(rand() % 256, rand() % 256, rand() % 256);
}

std::vector<cv::Scalar> RandomColors(int n) {
  std::vector<cv::Scalar> colors;
  colors.resize(n);
  for (int i = 0; i < n; ++i) {
    colors[i] = RandomColor();
  }
  return colors;
}

cv::Mat DrawImagesAndPoints(cv::InputArray image0,
                            cv::InputArray image1,
                            const std::vector<cv::Point2f> &points0,
                            const std::vector<cv::Point2f> &points1,
                            int point_radius,
                            int line_width) {
  CHECK(points0.size() == points1.size());
  auto colors = RandomColors(points0.size());
  return DrawImagesAndPoints(image0, image1,
                             points0, points1,
                             colors, colors,
                             point_radius, line_width);
}

// void WriteImageOrDie(const std::string &filename, const cv::Mat &image) {
//   LOG(INFO) << "Writing " << filename;
//   CHECK(cv::imwrite(filename, image))
//       << "Failed to write image \"" << filename << "\"";
// }

cv::Mat ReadImageOrDie(const std::string &filename) {
  auto image = cv::imread(filename);
  CHECK(!image.empty()) << "Failed to read image \"" << filename << "\"";
  return image;
}

} // furry
