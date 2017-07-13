#ifndef FURRY_SAC_ESTIMATOR_H
#define FURRY_SAC_ESTIMATOR_H

#include <vector>
#include <limits>
#include <iostream>
#include <Eigen/Dense>
#include <utility>

namespace furry
{

template <typename ModelDataType_,
          typename InputDataType_,
          typename ParamType_>
class SacModel
{
public:
  typedef InputDataType_ DataType;
  typedef ParamType_ ParamType;
  typedef ModelDataType_ ModelDataType;

  virtual void fit(const std::vector<DataType>& points) = 0;
  virtual double residual(const DataType& point) const = 0;
  ModelDataType data() const;
  void setParam(ParamType);
protected:
  ModelDataType data_;
  ParamType param_;
};

template <typename ModelDataType_,
          typename InputDataType_,
          typename ParamType_> std::ostream&
operator << (std::ostream& os,
             const SacModel<ModelDataType_, InputDataType_, ParamType_>& model);

class Line2DSacModel : public SacModel<Eigen::Vector3d, Eigen::Vector2d, int>
{
public:
  static const int kNumConstraints = 2;

  Line2DSacModel(ParamType param);
  virtual void fit(const std::vector<DataType>& points);
  virtual double residual(const DataType& point) const;
}; // Line2DSacModel

std::ostream& operator << (std::ostream& os, const Line2DSacModel& model);

class FundmentalMatrixSacModel
  : public SacModel<Eigen::Matrix3d,
                    std::pair<Eigen::Vector2d, Eigen::Vector2d>,
                    int>
{
public:
  static const int kNumConstraints = 8;

  FundmentalMatrixSacModel(ParamType param);
  virtual void fit(const std::vector<DataType>& points);
  virtual double residual(const DataType& point) const;
};

/**
 * Required member of SacModel:
 * DataType : data type of each element or constraint
 * ParamType : type of model parameters
 * kNumConstraints : number of data needed for fit
 * kdefaultParam : default model parameter
 * fit() : fit the model with required number of constraints
 * residual() : get the residual given the data
 */
template <typename SacModel>
class SacEstimator
{
public:
  typedef typename SacModel::DataType DataType;

  int numRounds() const;
  void setData(const std::vector<DataType>& data);
  SacModel model() const;
  virtual void run() = 0;
  virtual void reset();

protected:
  typename SacModel::ParamType model_param_ = typename SacModel::ParamType();
  SacModel model_ = SacModel(model_param_);
  std::vector<DataType> data_;
  int num_rounds_ = 0;
}; // SacEstimator

template <typename SacModel>
class RansacEstimator : public SacEstimator<SacModel>
{
public:
  //  typedef typename SacEstimator<SacModel>::DataType DataType;
  using typename SacEstimator<SacModel>::DataType;

  virtual void run();
  virtual void reset();
  void residualThreshold(double t);
  double residualThreshold() const;
  size_t numInliers() const;
  std::vector<bool> inlierMask() const;

protected:
  using SacEstimator<SacModel>::model_param_;
  using SacEstimator<SacModel>::model_;
  using SacEstimator<SacModel>::data_;
  using SacEstimator<SacModel>::num_rounds_;

  double residual_threshold_ = 1;
  size_t num_inliers_ = 0;
}; // RansacEstimator

} // furry

#include "sac_estimator_inl.hpp"

#endif // FURRY_SAC_ESTIMATOR_H
