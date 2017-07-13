const char *kUsageMessage = "Test segmentation methods and parameters";

#include <string>
#include <vector>
#include <cstdint>

#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>

#include "furry/3rdparty/slic/SLIC.h"
#include "furry/common/init.h"
#include "furry/common/clock.h"

DEFINE_string(in_image, "", "");
DEFINE_string(out_image, "", "");
DEFINE_int32(scale, 0, "");

// Mean Shift parameters
DEFINE_bool(use_ms, false, "");
DEFINE_double(ms_sp, 5, "");
DEFINE_double(ms_sr, 20, "");
DEFINE_int32(ms_max_level, 1, "");

// SLIC parameters
DEFINE_bool(use_slic, false, "");
DEFINE_int32(slic_k, 512, "");
DEFINE_int32(slic_m, 10, "");

using namespace cv;
using namespace std;

void DoMeanShift(const Mat &in_image, Mat *out_image) {
  pyrMeanShiftFiltering(in_image, *out_image,
                        FLAGS_ms_sp,
                        FLAGS_ms_sr,
                        FLAGS_ms_max_level);
}

void ReadSlicImage(const string &filename, unsigned int **pbuff, Size *size) {
  Mat image = imread(FLAGS_in_image);
  cvtColor(image, image, CV_BGR2RGB);
  *size = image.size();
  vector<Mat> channels(4);
  split(image, channels.data() + 1);
  channels[0] = Mat::ones(size->height, size->width, CV_8U) * 255;
  merge(channels, image);
  int num_ints = size->width * size->height;
  int num_bytes = num_ints * 4;
  *pbuff = new uint32_t[num_ints];
  memcpy(*pbuff, image.ptr(), num_bytes);
}

void Cv2Slic(const Mat &in_image, uint32_t **pbuff, Size *size) {
  Mat image;
  cvtColor(in_image, image, CV_BGR2RGB);
  *size = image.size();
  vector<Mat> channels(4);
  split(image, channels.data() + 1);
  channels[0] = Mat::ones(size->height, size->width, CV_8U) * 255;
  merge(channels, image);
  int num_ints = size->width * size->height;
  int num_bytes = num_ints * 4;
  *pbuff = new uint32_t[num_ints];
  memcpy(*pbuff, image.ptr(), num_bytes);
}

void SaveSlicImage(const string &filename,
                   const uint32_t *pbuff,
                   Size size) {
  Mat image(size, CV_8UC4);
  int num_bytes = size.width * size.height * 4;
  memcpy(image.ptr(), pbuff, num_bytes);
  vector<Mat> channels(4);
  split(image, channels.data());
  merge(channels.data() + 1, 3, image);
  cvtColor(image, image, CV_RGB2BGR);
  imwrite(filename, image);
}

void Slic2Cv(const uint32_t *pbuff, Size size, Mat *out_image) {
  out_image->create(size, CV_8UC4);
  int num_bytes = size.width * size.height * 4;
  memcpy(out_image->ptr(), pbuff, num_bytes);
  vector<Mat> channels(4);
  split(*out_image, channels.data());
  merge(channels.data() + 1, 3, *out_image);
  cvtColor(*out_image, *out_image, CV_RGB2BGR);
}

void DoSlic(const Mat &in_image, Mat *out_image) {
  // unsigned int (32 bits) to hold a pixel in ARGB format as follows:
  // from left to right,
  // the first 8 bits are for the alpha channel (and are ignored)
  // the next 8 bits are for the red channel
  // the next 8 bits are for the green channel
  // the last 8 bits are for the blue channel

  // ReadImage(pbuff, width, height);//YOUR own function to read an image into the ARGB format

  Size size;
  unsigned int* pbuff;
  // ReadSlicImage(FLAGS_in_image, &pbuff, &size);
  Cv2Slic(in_image, &pbuff, &size);
  int width = size.width;
  int height = size.height;

  //----------------------------------
  // Initialize parameters
  //----------------------------------
  int k = FLAGS_slic_k;//Desired number of superpixels.
  double m = FLAGS_slic_m;//Compactness factor. use a value ranging from 10 to 40 depending on your needs. Default is 10
  int* klabels = NULL;
  int numlabels(0);
  // string filename = FLAGS_out_image;
  // string savepath = "./";
  //----------------------------------
  // Perform SLIC on the image buffer
  //----------------------------------
  SLIC segment;
  segment.DoSuperpixelSegmentation_ForGivenNumberOfSuperpixels(pbuff, width, height, klabels, numlabels, k, m);
  // Alternately one can also use the function DoSuperpixelSegmentation_ForGivenStepSize() for a desired superpixel size
  //----------------------------------
  // Save the labels to a text file
  //----------------------------------
  // segment.SaveSuperpixelLabels(klabels, width, height, filename, savepath);
  //----------------------------------
  // Draw boundaries around segments
  //----------------------------------
  segment.DrawContoursAroundSegments(pbuff, klabels, width, height, 0xff0000);
  //----------------------------------
  // Save the image with segment boundaries.
  //----------------------------------
  // SaveSegmentedImageFile(pbuff, width, height);//YOUR own function to save an ARGB buffer as an image
  // SaveSlicImage(FLAGS_out_image, pbuff, size);
  Slic2Cv(pbuff, size, out_image);
  //----------------------------------
  // Clean up
  //----------------------------------
  if(pbuff) delete [] pbuff;
  if(klabels) delete [] klabels;
}

void ScaleImage(const Mat &src, int scale, Mat *dst) {
  *dst = src;
  for (int i = 0; i < scale; ++i) {
    pyrDown(*dst, *dst);
  }
}

int main(int argc, char *argv[]) {
  furry::Init(&argc, &argv, kUsageMessage);

  Mat in_image, out_image;
  in_image = imread(FLAGS_in_image);
  ScaleImage(in_image, FLAGS_scale, &in_image);

  furry::tic();
  if (FLAGS_use_ms) {
    DoMeanShift(in_image, &out_image);
  } else if (FLAGS_use_slic) {
    DoSlic(in_image, &out_image);
  } else {
    LOG(WARNING) << "Nothing is done";
    return 0;
  }
  furry::PrintEventTime("Segmentation");

  imwrite(FLAGS_out_image, out_image);

  return 0;
}
