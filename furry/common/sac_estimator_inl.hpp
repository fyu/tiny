#include "furry/common/rand.hpp"
#include "furry/common/debug.hpp"

namespace furry
{

////////////////////////////////////////////////////////////
// SacModel Inline
////////////////////////////////////////////////////////////

template <typename ModelDataType_,
          typename InputDataType_,
          typename ParamType_> inline ModelDataType_
SacModel<ModelDataType_, InputDataType_, ParamType_>::data() const
{
  return data_;
}

template <typename ModelDataType_,
          typename InputDataType_,
          typename ParamType_> inline void
SacModel<ModelDataType_, InputDataType_, ParamType_>::setParam(ParamType param)
{
  param_ = param;
}

template <typename ModelDataType_,
          typename InputDataType_,
          typename ParamType_> std::ostream&
operator << (std::ostream& os,
             const SacModel<ModelDataType_, InputDataType_, ParamType_>& model)
{
  os << model.data();
  return os;
}

////////////////////////////////////////////////////////////
// SacEstimator Inline
////////////////////////////////////////////////////////////

template <typename SacModel> inline void
SacEstimator<SacModel>::setData(const std::vector<DataType>& data)
{
  data_ = data;
  reset();
}

template <typename SacModel> inline int
SacEstimator<SacModel>::numRounds() const
{
  return num_rounds_;
}

template <typename SacModel> inline void
SacEstimator<SacModel>::reset()
{
  model_ = SacModel(model_param_);
  num_rounds_ = 0;
}

template <typename SacModel> inline SacModel
SacEstimator<SacModel>::model() const
{
  return model_;
}

////////////////////////////////////////////////////////////
// RansacEstimator Inline
////////////////////////////////////////////////////////////

template <typename SacModel> inline void
RansacEstimator<SacModel>::run()
{
  furry_debug_assert(data_.size() >= SacModel::kNumConstraints);
  //std::vector<DataType> samples = sample_unique(data_, SacModel::kNumConstraints);
  std::vector<int> sample_indices = deal(data_.size(), SacModel::kNumConstraints);
  std::vector<DataType> samples(SacModel::kNumConstraints);
  for (int i = 0; i < sample_indices.size(); ++i)
    samples[i] = data_[sample_indices[i]];
  SacModel new_model(model_param_);
  new_model.fit(samples);
  int new_num_inliers = 0;
  for (int i = 0; i < data_.size(); ++i)
  {
    if (new_model.residual(data_[i]) < residual_threshold_)
    {
      ++new_num_inliers;
    }
  }
  if (num_rounds_ == 0 || new_num_inliers >= num_inliers_)
  {
    model_ = new_model;
    num_inliers_ = new_num_inliers;
  }
  num_rounds_++;
}

template <typename SacModel> inline void
RansacEstimator<SacModel>::reset()
{
  SacEstimator<SacModel>::reset();
  num_inliers_ = 0;
}

template <typename SacModel> inline void
RansacEstimator<SacModel>::residualThreshold(double t)
{
  residual_threshold_ = t;
}

template <typename SacModel> inline double
RansacEstimator<SacModel>::residualThreshold() const
{
  return residual_threshold_;
}

template <typename SacModel> inline size_t
RansacEstimator<SacModel>::numInliers() const
{
  return num_inliers_;
}

template <typename SacModel> inline std::vector<bool>
RansacEstimator<SacModel>::inlierMask() const
{
  std::vector<bool> mask(data_.size(), false);
  for (int i = 0; i < data_.size(); ++i)
    if (model_.residual(data_[i]) < residual_threshold_)
      mask[i] = true;
  return mask;
}

} // furry
