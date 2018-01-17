#ifndef FURRY_COMMON_ALGEBRA
#define FURRY_COMMON_ALGEBRA

#include <Eigen/Dense>

namespace furry
{

/*
 * compute the matrix rank
 */
template <typename Derived> int
rank(const Eigen::MatrixBase<Derived>& m)
{
  Eigen::FullPivLU<typename Derived::PlainObject> lu(m);
  return lu.rank();
}

/*
 * find matrix m' that minimizes the Frobenius norm ||m - m'||
 * subject to the condition det m' = 0.
 */
template <typename Derived> Derived
closest_singular(const Eigen::MatrixBase<Derived>& m_)
{
  typedef Eigen::Matrix<typename Derived::Scalar,
                        Eigen::Dynamic,
                        Eigen::Dynamic>
      SVDInputType;
  SVDInputType m = m_;
  Eigen::JacobiSVD<SVDInputType> svd(m, Eigen::ComputeThinV | Eigen::ComputeThinU);
  auto singular_values = svd.singularValues();
  // std::cout << "sigular values: " << svd.singularValues().transpose()
  //           << " " << svd.singularValues()(0) / svd.singularValues()(2) << "\n\n";
  singular_values[singular_values.size()-1] = 0;
  return svd.matrixU() * singular_values.asDiagonal() * svd.matrixV().adjoint();
}

} // furry

#endif // FURRY_COMMON_ALGEBRA
