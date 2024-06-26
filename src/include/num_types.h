#pragma once
#ifndef __NUM_TYPES_H__
#define __NUM_TYPES_H__

#include<Eigen/Core>
#include<Eigen/Dense>

namespace types {
	using vec3 = Eigen::Vector3d;
	using vecd = Eigen::VectorXd;
	using mat3 = Eigen::Matrix3d;
	using mat2 = Eigen::Matrix2d;
	using matd = Eigen::MatrixXd;
}

#endif __NUM_TYPES_H__