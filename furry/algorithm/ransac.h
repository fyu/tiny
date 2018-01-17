#ifndef FURRY_ALGORITHM_RANSAC_H
#define FURRY_ALGORITHM_RANSAC_H

#include <functional>
#include <vector>

namespace furry {

template <typename Model, typename Data>
class Ransac {
 public:

  typedef std::function<void (const std::vector<Data> &data,
                              std::vector<Model> *models)> ModelFunc;
  typedef std::function<bool (const std::vector<Data> &data)> DegenFunc;
  typedef std::function<double (const Model &model, const Data &data)> DistFunc;

  struct Options {
    int num_samples;
    int max_num_trials = 2000;
    int max_num_data_trials = 100;
    float dist_threshold;
    float confidence = 0.99;
    ModelFunc model_func;
    DegenFunc degen_func;
    DistFunc dist_func;
  };

  Ransac(const Options &options) : options_(options) {}

  bool Fit(const std::vector<Data> &data,
           Model *best_model,
           std::vector<int> *best_inliers,
           int *final_num_trials = nullptr) const;

 private:

  void FindInliers(const Model &model,
                   const std::vector<Data> &data,
                   const DistFunc &dist_func,
                   float dist_threshold,
                   std::vector<int> *inliers) const;
  void FindBestModel(const std::vector<Model> &models,
                     const std::vector<Data> &data,
                     const DistFunc &dist_func,
                     float dist_threshold,
                     Model *best_model,
                     std::vector<int> *best_inliers) const;

 private:
  Options options_;

};

} // furry

#include "ransac-inl.h"

#endif
