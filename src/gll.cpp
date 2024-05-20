#include <gll.h>
#include <cassert>

namespace gll
{
	// derivative[j][point]
	shape_func::shape_func()
	{
		derivative.resize(max_p - 1);
		for (int order = 0; order < max_p - 1; order++)
		{
			derivative[order].resize(order + 2);
			for (int l_j = 0; l_j < order + 2; l_j++)
			{
				derivative[order][l_j].resize(order + 2);
				for (int point = 0; point < order + 2; point++)
				{
					for (int point_i = 0; point_i < order + 2; point_i++)
					{
						if (l_j == point_i)
						{
							continue;
						}
						double tmp = 1 / (points[order][l_j] - points[order][point_i]);

						for (int point_m = 0; point_m < order + 2; point_m++)
						{
							if (point_m == point_i || point_m == l_j)
							{
								continue;
							}
							tmp *= (points[order][point] - points[order][point_m]) / (points[order][l_j] - points[order][point_m]);
						}
						derivative[order][l_j][point] += tmp;
					}
				}
			}
		}
	}
	double shape_func::dl(int order, int shape_index, int point_index) {
		assert(order + 1< gll::max_p);
		return derivative[order + 1][shape_index][point_index];
	}
	double shape_func::l(int order, int shape_index, int point_index) {
		assert(order + 1< gll::max_p);
		return (shape_index == point_index) ? 1. : 0.;
	}
	shape_func func;
}