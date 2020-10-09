#ifndef _TROP_BASIC_HH_
#define _TROP_BASIC_HH_

#include <cmath>
#include <vector>
#include <Eigen/Core>
#include <Eigen/SparseCore>

using std::fabs;
using std::round;
using std::vector;
using std::string;

typedef Eigen::Matrix<double, Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> RowMatrix;
typedef Eigen::Matrix<double, Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor> ColMatrix;

typedef Eigen::Matrix<double, 1,Eigen::Dynamic> RowVector;
typedef Eigen::Matrix<double, Eigen::Dynamic,1> ColVector;

typedef Eigen::Matrix<int, 1,Eigen::Dynamic> VectorI;

typedef Eigen::SparseMatrix<double,Eigen::RowMajor> SparseRowMatrix;

typedef long long VolType;

#endif
