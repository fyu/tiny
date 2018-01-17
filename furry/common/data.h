#ifndef FURRY_COMMON_DATA_H
#define FURRY_COMMON_DATA_H

#include <algorithm>
#include <functional>
#include <iostream>
#include <vector>
#include <numeric>
#include <string>
#include <cmath>

namespace furry {

class DataStat {
 public:
  DataStat();
  DataStat(const std::vector<double> &data,
           std::function<bool (double, double)> comp_func = std::less<double>());
  DataStat(DataStat &&ds);

  void SetTitle(const std::string &title);
  void SetTitle(std::string &&title);
  void SetSimpleOutput(bool simple);

  DataStat& operator = (DataStat &&ds);

 protected:
  std::vector<double> data_;
  std::function<bool (double, double)> comp_func_;
  double mean_;
  double sigma_;
  std::string title_;
  bool simple_output_ = false;

  friend std::ostream& operator << (std::ostream &os, const DataStat &ds);
};

std::ostream& operator << (std::ostream &os, const DataStat &ds);

////////////////////////////////////////////////////////////
// DataStat Inline
////////////////////////////////////////////////////////////

inline DataStat::DataStat() {}

inline DataStat::DataStat(const std::vector<double> &data,
                          std::function<bool (double, double)> comp_func)
    : data_(data), comp_func_(comp_func) {
  if (data.size() != 0) {
    std::sort(data_.begin(), data_.end(), comp_func_);
    mean_ = std::accumulate(data_.begin(), data_.end(), 0.0) / data_.size();
    sigma_ = sqrt(std::accumulate(data_.begin(), data_.end(), 0.0,
                                  [&](double a, double b) {
                                    return a + (b - mean_) * (b - mean_);
                                  }) / data_.size());
  }
}

inline DataStat::DataStat(DataStat &&ds) {
  *this = std::move(ds);
}

inline void DataStat::SetTitle(const std::string &title) {
  title_ = title;
}

inline void DataStat::SetTitle(std::string &&title) {
  title_ = std::move(title);
}

inline void DataStat::SetSimpleOutput(bool simple) {
  simple_output_ = simple;
}

inline DataStat& DataStat::operator = (DataStat &&ds) {
  if (this != &ds) {
    data_ = std::move(ds.data_);
    title_ = std::move(ds.title_);
    mean_ = ds.mean_;
    sigma_ = ds.sigma_;
    comp_func_ = ds.comp_func_;
  }
  return *this;
}

double GetMean(const std::vector<double> &v);
double GetStd(const std::vector<double> &v);
double Max(const std::vector<double> &v);

} // furry

#endif
