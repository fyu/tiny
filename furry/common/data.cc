#include "furry/common/data.h"

#include <iomanip>
#include <cmath>
#include <limits>

namespace furry {

std::ostream& operator << (std::ostream &os, const DataStat &ds) {
  if (ds.data_.size() == 0) {
    os << "[DataStat] No data\n";
    return os;
  }
  if (ds.title_.size() > 0) {
    os << "\n============================================================\n"
       << "\t" << ds.title_ << '\n'
       << "============================================================\n\n";
  }
  os << "# Items: " << ds.data_.size() << '\n';
  if (ds.simple_output_) {
    os << ds.mean_ << '\n'
       << ds.sigma_ << '\n';
    for (int i = 0; i < 100; i += 10)
      os << ds.data_[ds.data_.size() * i / 100.0] << '\n';
    os << ds.data_[ds.data_.size() - 1] << '\n';
  } else {
    os << "Mean: " << ds.mean_ << '\n'
       << "Sigma: " << ds.sigma_ << '\n';
    for (int i = 0; i < 100; i += 10) {
      os << std::setw(3) << std::setfill(' ') << 100 - i
         << std::setw(0) << "% percentile: "
         << ds.data_[ds.data_.size() * i / 100.0]
         << '\n';
    }
    os << "  0% percentile: " << ds.data_[ds.data_.size() - 1] << '\n';
  }
  return os;
}

double GetMean(const std::vector<double> &v) {
  if (v.size() == 0) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  double mean = 0;
  for (double a : v) {
    mean += a;
  }
  return mean / v.size();
}

double GetStd(const std::vector<double> &v) {
  double mean = GetMean(v);
  double s = 0;
  for (double a : v) {
    s += (a - mean) * (a - mean);
  }
  return sqrt(s / (v.size() - 1));
}

double Max(const std::vector<double> &v) {
  if (v.size() == 0) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  double m = v[0];
  for (double a : v) {
    if (a > m) m = a;
  }
  return m;
}

}
