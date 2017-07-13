#include "plane_sweeping.h"

#include "furry/common/cv.h"
#include "furry/common/memory.h"

namespace furry {
namespace tiny {

bool PlaneSweepingData::WriteBin(FILE *file) const {
  fwrite(&min_depth, sizeof(double), 1, file);
  fwrite(&max_depth, sizeof(double), 1, file);
  fwrite(&patch_radius, sizeof(int), 1, file);
  fwrite(&scale, sizeof(int), 1, file);
  fwrite(&num_samples, sizeof(int), 1, file);
  fwrite(intrinsic_matrix.data(), sizeof(double), 9, file);
  return WriteImage(image, file);
}

bool PlaneSweepingData::ReadBin(FILE *file) {
  fread(&min_depth, sizeof(double), 1, file);
  fread(&max_depth, sizeof(double), 1, file);
  fread(&patch_radius, sizeof(int), 1, file);
  fread(&scale, sizeof(int), 1, file);
  fread(&num_samples, sizeof(int), 1, file);
  fread(intrinsic_matrix.data(), sizeof(double), 9, file);
  return ReadImage(file, &image);
}

bool PlaneSweepingData::WriteTxt(const std::string &filename) const {
  auto file = ScopedPtr<FILE>(fopen(filename.c_str(), "w"));
  fprintf(file, "%12s : %-lf\n", "min_depth", min_depth);
  fprintf(file, "%12s : %-lf\n", "max_depth", max_depth);
  fprintf(file, "%12s : %-d\n", "patch_size", patch_radius * 2 + 1);
  fprintf(file, "%12s : %-d\n", "scale", scale);
  fprintf(file, "%12s : %-d\n", "num_planes", num_samples);
  return true;
}

} // tiny
} // furry
