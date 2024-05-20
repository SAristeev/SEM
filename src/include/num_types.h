#pragma once
#ifndef __NUM_TYPES_H__
#define __NUM_TYPES_H__

#include<Eigen/Core>

namespace types {
	using vec3 = Eigen::Vector3d;
	using mat3 = Eigen::Matrix3d;
	using matd = Eigen::MatrixXd;
}

#endif __NUM_TYPES_H__