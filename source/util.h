// util.h

#pragma once

#include <Eigen/Dense>
#include <complex>
#include <string>

typedef std::complex<double> Complex;
typedef Eigen::Vector3<double> Vec3;
typedef Eigen::Vector4<double> Vec4;

/// @brief Same as sprintf but generates a std::string
/// @tparam ...Args Template for string arguments
/// @param format Format string (same as sprintf)
/// @param ...args String arguments
/// @return A formatted std::string
template <typename... Args>
std::string string_format(const char* format, Args... args) {
  int size = std::snprintf(nullptr, 0, format, args...) + 1;
  assert(size > 0);
  char buf[size];
  sprintf(buf, format, args...);
  return std::string(buf);
}

/// @brief Check if 2 numbers are almost equal
/// @param x1 First number
/// @param x2 Second number
/// @param eps Numbers are almost equal if their difference is within +/-
/// epsilon
/// @return true if almost equal, false otherwise
bool AlmostEq(const double x1, const double x2, double eps = 1.0e-14) {
  return fabs(x1 - x2) < eps;
}

/// @brief Check if 2 vectors are almost equal
/// @param v1 First vector
/// @param v2 Second vector
/// @param eps Vectors are almost equal if v1.v2 is within 1 +/- epsilon
/// @return true if almost equal, false otherwise
bool AlmostEq(const Eigen::VectorXd& v1, const Eigen::VectorXd& v2,
              double eps = 1.0e-14) {
  return fabs(v1.dot(v2) - 1.0) < eps;
}

/// @brief Chop a floating point number that is close to zero
/// @param x Number to chop
/// @param eps Numbers with absolute values less than epsilon will be chopped
/// @return Chopped number
double Chop(double x, double eps = 1.0e-14) {
  if (fabs(x) < eps) return 0.0;
  return x;
}

/// @brief Convert a 3-vector to a hashable string
/// @param v 3-Vector
/// @param n_digits Number of digits to trim each coordinate value to
/// @return Hashable string that can be used to identify this vector
std::string Vec3ToString(const Vec3 v, int n_digits = 10) {
  double eps = pow(10.0, double(-n_digits));
  double v_x = Chop(v.x(), eps);
  double v_y = Chop(v.y(), eps);
  double v_z = Chop(v.z(), eps);
  return string_format("(%+.*f,%+.*f,%+.*f)", n_digits, v_x, n_digits, v_y,
                       n_digits, v_z);
}

/// @brief Convert a 4-vector to a hashable string
/// @param v 4-Vector
/// @param n_digits Number of digits to trim each coordinate value to
/// @return Hashable string that can be used to identify this vector
std::string Vec4ToString(const Vec4 v, int n_digits = 10) {
  double eps = pow(10.0, double(-n_digits));
  double v_x = Chop(v.x(), eps);
  double v_y = Chop(v.y(), eps);
  double v_z = Chop(v.z(), eps);
  double v_w = Chop(v.w(), eps);
  return string_format("(%+.*f,%+.*f,%+.*f,%+.*f)", n_digits, v_x, n_digits,
                       v_y, n_digits, v_z, n_digits, v_w);
}
