#pragma once
#ifndef __PARAMETRIC_HEX_H__
#define __PARAMETRIC_HEX_H__

#include "num_types.h"
#include <array>
#include <vector>
namespace pre
{
	//                                              0 1 2 3 4 5
	//                                              ! ! ! ! ! !
	static std::array<int, 6> faces_fidesys2vtk = { 4,0,3,1,2,5 };

	struct vtk_faces
	{
		std::array<std::array<int, 4>, 6> faces;
		vtk_faces()
		{
			{
				faces[0] = { 0, 1, 4, 5 };
				faces[1] = { 3, 2, 7, 6 };

				faces[2] = { 0, 3, 4, 7 };
				faces[3] = { 1, 2, 5, 6 };
				
				faces[4] = { 0, 1, 3, 2 };
				faces[5] = { 4, 5, 7, 6 };
			}
		}
	};

	using namespace types;
	struct boundary_indexes
	{
		std::array<std::array<int, 2>, 12> edges;
		std::array<std::array<int, 4>, 6> faces;
		boundary_indexes()
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
				//incorrect
				faces[0] = { 0, 3, 4, 7 };
				faces[1] = { 1, 2, 5, 6 };
				faces[2] = { 0, 1, 4, 5 };
				faces[3] = { 3, 2, 7, 6 };
				faces[4] = { 0, 1, 3, 2 };
				faces[5] = { 4, 5, 7, 6 };
			}
		}
	};
	
	vec3 get_outer_normal(int face_id);
	std::array<int, 3> get_local_index(int order, int i);
	void get_face(int order, int id, std::vector<int>& face_id);
}
#endif