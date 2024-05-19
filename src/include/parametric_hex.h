#pragma once
#ifndef __PARAMETRIC_HEX_H__
#define __PARAMETRIC_HEX_H__
namespace pre 
{
	struct local_indexing
	{
		std::array<std::array<int, 2>, 12> edges;
		std::array<std::array<int, 4>, 6> faces;
		local_indexing() 
		{
			{
				edges[0] = { 0, 1 };
				edges[1] = { 1, 2 };
				edges[2] = { 3, 2 };
				edges[3] = { 0, 3 };

				edges[4] = { 4, 5 };
				edges[5] = { 5, 6 };
				edges[6] = { 7, 6 };
				edges[7] = { 4, 7 };

				edges[8] = { 0, 4 };
				edges[9] = { 1, 5 };
				edges[10] = { 2, 6 };
				edges[11] = { 3, 7 };
			}
			{
				faces[0] = { 0, 3, 4, 7 };
				faces[1] = { 1, 2, 5, 6 };
				faces[2] = { 0, 1, 4, 5 };
				faces[3] = { 3, 2, 7, 6 };
				faces[4] = { 0, 1, 3, 2 };
				faces[5] = { 4, 5, 7, 6 };
			}
		}
		std::tuple<int, int, int> get_indx(int id, int order) 
		{

		}
	};
}
#endif