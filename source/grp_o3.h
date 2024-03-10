// grp_o3.h

#pragma once

#include <Eigen/Dense>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <map>
#include <string>
#include <vector>

typedef Eigen::Quaternion<double> Quat;
typedef Eigen::Vector3<double> Vec3;

class GrpElemO3 {
 public:
  GrpElemO3();
  GrpElemO3(Quat q, int inv);
  GrpElemO3(double w, double x, double y, double z, int inv);
  Vec3 operator*(Vec3 const& v) const;
  GrpElemO3 operator*(GrpElemO3 const& g) const;
  void ReadGrpElem(FILE* file);

  Quat q;   // SO(3) rotation in quaternion representation
  int inv;  // -1 if this group element includes an inversion
};

GrpElemO3::GrpElemO3() {
  this->q = Quat(0.0, 0.0, 0.0, 1.0);
  this->inv = 1;
}

/// @brief O(3) group element constructor
/// @param q SO(3) rotation in quaternion representation
/// @param inv true if this group element includes an inversion
GrpElemO3::GrpElemO3(Quat q, int inv) {
  this->q = q;
  this->inv = inv;
}

/// @brief O(3) group element constructor
/// @param w Quaternion real part
/// @param x Quaternion "i" coefficient
/// @param y Quaternion "j" coefficient
/// @param z Quaternion "k" coefficient
/// @param inv true if this group element includes an inversion
GrpElemO3::GrpElemO3(double w, double x, double y, double z, int inv) {
  q = Quat(w, x, y, z);
  this->inv = inv;
}

/// @brief Rotate a vector by this group element
/// @param v Vector to rotate
/// @return Rotated vector
Vec3 GrpElemO3::operator*(Vec3 const& v) const { return q * v * double(inv); }

/// @brief Multiply by another group element
/// @param g Group element to multiply
/// @return Product of this element on the left times @p g on the right
GrpElemO3 GrpElemO3::operator*(GrpElemO3 const& g) const {
  return GrpElemO3(q * g.q, inv * g.inv);
}

void GrpElemO3::ReadGrpElem(FILE* file) {
  int index;
  double q_w, q_x, q_y, q_z;
  fscanf(file, "%d %d %lf %lf %lf %lf\n", &index, &inv, &q_w, &q_x, &q_y, &q_z);
  q = Quat(q_w, q_x, q_y, q_z);
}