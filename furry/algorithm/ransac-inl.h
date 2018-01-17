#include <algorithm>
#include <cmath>
#include <utility>

#include <glog/logging.h>

#include "furry/common/rand.h"
#include "furry/common/sequence.h"

namespace furry {

template <typename Model, typename Data> inline
bool Ransac<Model, Data>::Fit(const std::vector<Data> &data,
                              Model *best_model,
                              std::vector<int> *best_inliers,
                              int *final_num_trials) const {
  if (data.size() < options_.num_samples)
    return false;
  int num_trials = 0;
  int total_trials = 1;
  int total_count = 1;
  int num_points = data.size();
  std::vector<int> sample_indexes(options_.num_samples);
  std::vector<Model> models;
  Model model;
  std::vector<int> inliers;
  double eps = 1e-15;

  best_inliers->resize(0);

  while (total_trials > num_trials) {
    bool degenerate = true;
    int count = 1;
    while (degenerate) {
      SampleUnique(num_points, options_.num_samples, sample_indexes.data());
      auto samples = Slice(data, sample_indexes);
      degenerate = options_.degen_func(samples);
      if (!degenerate) {
        options_.model_func(samples, &models);
        if (models.size() == 0)
          degenerate = true;
      }
      count = count + 1;
      if (count > options_.max_num_data_trials) {
        LOG(ERROR) << "Unable to select a nondegenerate data set";
        return false;
      }
    }
    FindBestModel(models, data, options_.dist_func, options_.dist_threshold,
                  &model, &inliers);
    if (inliers.size() > best_inliers->size()) {
      swap(inliers, *best_inliers);
      *best_model = model;

      double frac_inliers = (double)best_inliers->size() / num_points;
      double p_no_outliers = 1 - pow(frac_inliers, options_.num_samples);
      // avoid division by -inf or 0
      p_no_outliers = fmin(fmax(eps, p_no_outliers), 1 - eps);
      total_trials = log(1 - options_.confidence) / log(p_no_outliers);
    }

    ++total_count;

    if (total_count > options_.max_num_trials) {
      // LOG(WARNING) << "Ransac reached max number of "
      //              << options_.max_num_trials << " trials";
      break;
    }
  }

  if (final_num_trials) {
    *final_num_trials = total_count;
  }

  return true;
}

template <typename Model, typename Data> inline
void Ransac<Model, Data>::FindInliers(const Model &model,
                                      const std::vector<Data> &data,
                                      const DistFunc &dist_func,
                                      float dist_threshold,
                                      std::vector<int> *inliers) const {
  inliers->resize(0);
  for (size_t i = 0; i < data.size(); ++i) {
    if (dist_func(model, data[i]) < dist_threshold) {
      inliers->push_back(i);
    }
  }
}

template <typename Model, typename Data> inline
void Ransac<Model, Data>::FindBestModel(const std::vector<Model> &models,
                                        const std::vector<Data> &data,
                                        const DistFunc &dist_func,
                                        float dist_threshold,
                                        Model *best_model,
                                        std::vector<int> *best_inliers) const {
  std::vector<int> inliers;
  for (size_t i = 0; i < models.size(); ++i) {
    FindInliers(models[i], data, dist_func, dist_threshold, &inliers);
    if (inliers.size() > best_inliers->size()) {
      std::swap(inliers, *best_inliers);
      *best_model = models[i];
    }
  }
}


} // furry
