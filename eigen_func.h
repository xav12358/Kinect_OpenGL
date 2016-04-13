#ifndef EIGEN_FUNC_H
#define EIGEN_FUNC_H


#include <include/Eigen/Dense>

Eigen::Vector3f crossE(Eigen::Vector3f u, Eigen::Vector3f v);
Eigen::Affine3d create_rotation_matrix(double ax, double ay, double az);
Eigen::Matrix4f exp(const Eigen::VectorXf& mu);
void            rodrigues_so3_exp(const Eigen::Vector3f w, const float A, const float B, Eigen::Matrix3f &R);

#endif // EIGEN_FUNC_H
