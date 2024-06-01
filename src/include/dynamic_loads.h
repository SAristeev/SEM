#pragma once
#ifndef __DYNAMIC_LOADS_H__
#define __DYNAMIC_LOADS_H__

#include <cmath>
#include <numbers>

inline double delta_load(double t, double t0, double A)
{
	return (std::abs(t - t0) < 1e-16) ? A : 0;	
}

inline double exp_load(double t, double A, double alpha)
{
	return A * std::exp(-1 * alpha * t);
}

inline double ricker_load(double t, double t0, double A, double w)
{
	return A * 
		   (1 - 2 * std::numbers::pi * std::numbers::pi * w * w * (t - t0) * (t - t0)) *
		   std::exp(-1 * std::numbers::pi * std::numbers::pi * w * w * (t - t0) * (t - t0));
}

inline double berlage_load(double t, double A, double w)
{
	double w0 = 2 * std::numbers::pi * w;
	double w1 = w0 / std::sqrt(3.);

	return
		A * w1 * w1 * std::exp(-1 * w1 * t) *
		(
			std::sin(w0 * t) *
				(-1 * t * t / w1 + t / (w1 * w1) + 1 / (w1 * w1 * w1))
			-
			std::cos(w0 * t) *
				std::sqrt(3.) * (t * t / w1 + t / (w1 * w1))
		) / 4;

	/*return
		A * w1 * w1 * std::exp(-1 * w1 * t) *
		(
			std::sin(w0 * t) *
			(-1 * t * t + t + 1) / (w1 * w1)
			-
			std::cos(w0 * t) *
			std::sqrt(3.) * (t * t + t ) / (w1 * w1)
			) / 4;*/
}


#endif //__DYNAMIC_LOADS_H__