#include "parametric_hex.h"

#include "num_types.h"
#include "gll.h"

#include <cassert>
#include <vector>
#include <unordered_map>

namespace pre
{
	using namespace types;

	static boundary_indexes bound;

	struct indexes_map
	{
		std::vector<std::unordered_map<int, std::array<int, 3>>> container;
		indexes_map()
		{
			std::vector<vec3> param_points(9);
			param_points[0] = { 0,0,0 };
			param_points[1] = { 1,0,0 };
			param_points[2] = { 1,1,0 };
			param_points[3] = { 0,1,0 };
			param_points[4] = { 0,0,1 };
			param_points[5] = { 1,0,1 };
			param_points[6] = { 1,1,1 };
			param_points[7] = { 0,1,1 };
			param_points[8] = { 0.5,0.5,0.5 };


			container.resize(gll::max_p - 1);
			
			for (int order = 0; order < gll::max_p - 1; order++)
			{
				int shift = 0;
				int nodes = order + 2;
				// vertices
				{
					int p = nodes - 1;
					container[order][0] = { 0,0,0 };
					container[order][1] = { p,0,0 };
					container[order][2] = { p,p,0 };
					container[order][3] = { 0,p,0 };
					container[order][4] = { 0,0,p };
					container[order][5] = { p,0,p };
					container[order][6] = { p,p,p };
					container[order][7] = { 0,p,p };
					shift = 8;
				}
				// edges
				{
					auto add_edge_loc = [&](int x_min_idx, int x_max_idx)
						{
							for (int i = 1; i < nodes - 1; i++)
							{							
								vec3 param_vec = param_points[x_max_idx] - param_points[x_min_idx];
								std::array<int, 3> tmp = { 0,0,0 };
								int forward_indx = 0;
								for (int dim = 0; dim < 3; dim++)
								{
									if (abs(param_vec[dim]) > 1e-16) { forward_indx = dim; break; }
								}
								bool is_right = (param_points[x_max_idx] - param_points[x_min_idx])[forward_indx] > 0;
								tmp[forward_indx] = is_right ? i : (nodes - 1 - i);

								vec3 orientation = param_points[x_min_idx] - param_points[8];
								tmp[(forward_indx + 1) % 3] = orientation[(forward_indx + 1) % 3] > 0 ? (nodes - 1) : 0;
								tmp[(forward_indx + 2) % 3] = orientation[(forward_indx + 2) % 3] > 0 ? (nodes - 1) : 0;
								container[order][shift] = { tmp[0], tmp[1], tmp[2] };
								shift++;
							}
						};

					for (auto& e : bound.edges)
					{
						add_edge_loc(e[0], e[1]);
					}
				}
				// faces
				{
					//  y
					//	^	 lu    ru
					//	|	 2. . .3
					//	|	   .   
					//	|	     . 
					//	|	 0. . .1
					//  |    ld     rd
					//  +-------------->x
					auto add_face_loc = [&](int ld_indx, int rd_indx, int lu_indx, int ru_indx)
						{
							for (int y_ind = 1; y_ind < nodes - 1; y_ind++)
							{
								for (int x_ind = 1; x_ind < nodes - 1; x_ind++)
								{
									// bug in this function
									vec3 param_vec_x = param_points[rd_indx] - param_points[ld_indx];
									vec3 param_vec_y = param_points[lu_indx] - param_points[ld_indx];
									std::array<int, 3> tmp = { 0,0,0 };
									int forward_indx_x = 0;
									int forward_indx_y = 0;
									for (int dim = 0; dim < 3; dim++)
									{
										if (abs(param_vec_x[dim]) > 1e-16) { forward_indx_x = dim; break; }
									}
									for (int dim = 0; dim < 3; dim++)
									{
										if (abs(param_vec_y[dim]) > 1e-16) { forward_indx_y = dim; break; }
									}

									bool is_right_x = (param_points[rd_indx] - param_points[ld_indx])[forward_indx_x] > 0;
									bool is_right_y = (param_points[lu_indx] - param_points[ld_indx])[forward_indx_y] > 0;

									tmp[forward_indx_x] = is_right_x ? x_ind : (nodes - 1 - x_ind);
									tmp[forward_indx_y] = is_right_y ? y_ind : (nodes - 1 - y_ind);

									vec3 orientation = param_points[ld_indx] - param_points[8];
									// bug
									tmp[3 - forward_indx_x - forward_indx_y] = orientation[3 - forward_indx_x - forward_indx_y] > 0 ? (nodes - 1) : 0;

									container[order][shift] = { tmp[0], tmp[1], tmp[2] };
									shift++;
								}
							}
						};

					for (auto& f : bound.faces)
					{
						add_face_loc(f[0], f[1], f[2], f[3]);
					}

				}
				// volume
				{
					for (int z_ind = 1; z_ind < nodes - 1; z_ind++)
					{
						for (int y_ind = 1; y_ind < nodes - 1; y_ind++)
						{
							for (int x_ind = 1; x_ind < nodes - 1; x_ind++)
							{
								std::array<int, 3> tmp = { x_ind, y_ind, z_ind };
								container[order][shift] = { tmp[0], tmp[1], tmp[2] };
								shift++;
							}
						}
					}
				}
			}
		}
	};

	static indexes_map indx_map;

	// \xi \eta \zeta
	std::array<int, 3> get_local_index(int order, int i)
	{
		return indx_map.container[order - 1][i];
	}

	vec3 get_outer_normal(int face_id)
	{
		switch (face_id)
		{
		case 0:
			return { 0, -1, 0 };
		case 1:
			return { 0, 1, 0 };
		case 2:
			return { -1, 0, 0 };
		case 3:
			return { 1, 0, 0 };
		case 4:
			return { 0, 0, -1 };
		case 5:	  
			return { 0, 0, 1 };
		default:
			assert(0);
			break;
		}
	}

	void get_face(int order, int face_id, std::vector<int>& face_idxs) 
	{
		int nodes3 = (order + 1) * (order + 1) * (order + 1);
		face_idxs.resize((order + 1) * (order + 1));
		std::fill(face_idxs.begin(), face_idxs.end(), 0);
		int id = 0;
		for (int i = 0; i < nodes3; ++i) 
		{
			std::array<int, 3> loc = get_local_index(order, i);
			switch (face_id)
			{
			case 0:
				if (loc[1] == 0){
					face_idxs[id] = i;
					id++;
				}
				break;
			case 1:
				if (loc[1] == order) {
					face_idxs[id] = i;
					id++;
				}
				break;
			case 2:
				if (loc[0] == 0) {
					face_idxs[id] = i;
					id++;
				}
				break;
			case 3:
				if (loc[0] == order) {
					face_idxs[id] = i;
					id++;
				}
				break;
			case 4:
				if (loc[2] == 0) {
					face_idxs[id] = i;
					id++;
				}
				break;
			case 5:
				if (loc[2] == order) {
					face_idxs[id] = i;
					id++;
				}
				break;
			default:
				break;
			}
		}
	}

	void get_nodeset(int order, std::array<int, 4> nodeset, std::vector<int>& face_idxs) 
	{
		for (int face = 0; face < 6; face++) 
		{
			std::array<bool, 4> find = { false, false, false, false };
			for (int point = 0; point < 4; point++) 
			{
				if (std::find(bound.faces[face].begin(), bound.faces[face].end(), nodeset[point]) != bound.faces[face].end()) 
				{
					find[point] = true;
				}
			}
			if (find[0] && find[1] && find[2] && find[3]) 
			{
				get_face(order, face, face_idxs);
				break;
			}
		}
	}
}