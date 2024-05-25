#include "fc.h"
#include "gll.h"
#include "parametric_hex.h"

#include <fstream>
#include <algorithm>	

#include <base64.h>
#include <json/json.hpp>

namespace pre {
	static boundary_indexes bound;

	void BC::from_nodeset(const UnstructedMesh& original_mesh, const UnstructedMesh& computational_mesh)
	{

		if (original_mesh.spectral_elems) {
			std::array<std::array<int, 8>, 6> face_mask;
			{
				face_mask[0] = { 0, 3, 4, 7, 11, 16, 15, 19 };
				face_mask[1] = { 1, 2, 5, 6, 9, 13, 15, 17 };
				face_mask[2] = { 0, 1, 4, 5, 8, 12, 16, 17 };
				face_mask[3] = { 2, 3, 6, 7, 10, 14, 18, 19 };
				face_mask[4] = { 0, 1, 2, 3, 8, 9, 10, 11};
				face_mask[5] = { 4, 6, 7, 8, 12, 13, 14, 15};
			}

			// elem, loc_id
			std::map<int, std::vector<int>> apply_to_map;
			for (int i = 0; i < apply_to.size(); i++) 
			{
				apply_to[i] = original_mesh.map_node_numeration.at(apply_to[i]);
			}

			for (int elem_id = 0; elem_id < original_mesh.elemids.size(); elem_id++) 
			{
				for (int point = original_mesh.elem_shifts[elem_id]; point < original_mesh.elem_shifts[elem_id + 1]; point++) 
				{
					for (int i = 0; i < apply_to.size(); i++) 
					{
						auto find = std::find(apply_to.begin(), apply_to.end(), original_mesh.elems[point]);
						if (find != apply_to.end()) 
						{
							auto unique = std::find(apply_to_map[elem_id].begin(), apply_to_map[elem_id].end(), original_mesh.elems[point]);
							if (unique == apply_to_map[elem_id].end())
							{
								apply_to_map[elem_id].push_back(*find);
							}
						}
					}
				}
			}
			std::vector<int> new_apply_to;
			for (int elem_id = 0; elem_id < original_mesh.elemids.size(); elem_id++) 
			{
				if(apply_to_map[elem_id].size() >= 8) 
				{
					for (int face = 0; face < 6; face++) 
					{
						std::array<bool, 8> finded = {false, false, false, false, false, false, false, false };
						for(int point = 0; point < 8; point++)
						{
							finded[point] = std::find(apply_to_map[elem_id].begin(), apply_to_map[elem_id].end(), original_mesh.elems[original_mesh.elem_shifts[elem_id] + face_mask[face][point]]) == apply_to_map[elem_id].end();
						}
						if (std::find(finded.begin(), finded.end(), false) == finded.end()) 
						{
							std::vector<int> tmp;
							get_face(original_mesh.order[elem_id], face, tmp);
							for (int i = 0; i < tmp.size(); i++) 
							{
								int to_add = computational_mesh.elems[computational_mesh.elem_shifts[elem_id] + tmp[i]];
								auto unique = std::find(new_apply_to.begin(), new_apply_to.end(), to_add);
								if (unique == new_apply_to.end())
								{
									new_apply_to.push_back(to_add);
								}
							}
						}
					}
				}
			}
			apply_to = new_apply_to;
		}
	}
	void BC::node_mapping(const UnstructedMesh& original_mesh, const UnstructedMesh& computational_mesh) 
	{
		for (int i = 0; i < apply_to.size(); i++)
		{
			vec3 orig_point = original_mesh.nodes[original_mesh.map_node_numeration.at(apply_to[i])];
			auto find = std::find_if(computational_mesh.nodes.begin(), computational_mesh.nodes.end(), [&](vec3 b) { return (orig_point - b).norm() < 1e-16; });
			if (find == computational_mesh.nodes.end())
			{
				assert(0);
			}
			else
			{
				apply_to[i] = find - computational_mesh.nodes.begin();
			}
		}

	}
	void fc::add_spectral_elem(int elem_id)
	{
		int original_offset = original_mesh.elem_shifts[elem_id];
		int order = computational_mesh.order[elem_id];
		int nodes = order + 1;
		int nodes3 = nodes * nodes * nodes;

		// vertices
		{
			for (int i = original_offset; i < original_offset + 8; i++)
			{
				int node = original_mesh.elems[i];
				vec3 point = original_mesh.nodes[node];
				auto find = std::find_if(computational_mesh.nodes.begin(), computational_mesh.nodes.end(), [&](vec3 b) { return (point - b).norm() < 1e-16; });
				if (find == computational_mesh.nodes.end())
				{
					computational_mesh.elems.push_back(computational_mesh.nodes.size());
					computational_mesh.nodes.push_back(point);
				}
				else
				{
					computational_mesh.elems.push_back(find - computational_mesh.nodes.begin());
				}
			}
		}
		// edges
		{
			auto add_edge = [&](int x_min_idx, int x_max_idx)
				{
					vec3 x_min = original_mesh.nodes[original_mesh.elems[x_min_idx + original_offset]];
					vec3 x_max = original_mesh.nodes[original_mesh.elems[x_max_idx + original_offset]];
					vec3 x_mid = (x_max + x_min) / 2;
					vec3 x_vec = x_max - x_mid;

					for (int i = 1; i < nodes - 1; i++)
					{
						vec3 point = gll::points[order - 1][i] * x_vec + x_mid;
						auto find = std::find_if(computational_mesh.nodes.begin(), computational_mesh.nodes.end(), [&](vec3 b) { return (point - b).norm() < 1e-16; });
						if (find == computational_mesh.nodes.end())
						{
							computational_mesh.elems.push_back(computational_mesh.nodes.size());
							computational_mesh.nodes.push_back(point);
						}
						else
						{
							computational_mesh.elems.push_back(find - computational_mesh.nodes.begin());
						}
					}
				};
			for (auto& e : bound.edges)
			{
				add_edge(e[0], e[1]);
			}
		}
		// faces
		{	
			auto add_face = [&](int ld_indx, int rd_indx, int lu_indx, int ru_indx)
				{
					//  y
					//	^	 lu    ru
					//	|	 2. . .3
					//	|	   .   
					//	|	     . 
					//	|	 0. . .1
					//  |    ld     rd
					//  +-------------->x
					
					vec3 ld = original_mesh.nodes[original_mesh.elems[ld_indx + original_offset]];
					vec3 rd = original_mesh.nodes[original_mesh.elems[rd_indx + original_offset]];
					vec3 lu = original_mesh.nodes[original_mesh.elems[lu_indx + original_offset]];
					vec3 ru = original_mesh.nodes[original_mesh.elems[ru_indx + original_offset]];

					vec3 d_mid = (rd + ld) / 2;
					vec3 l_mid = (lu + ld) / 2;

					vec3 vec_x = rd - d_mid;
					vec3 vec_y = lu - l_mid;

					vec3 mid = (lu + ld + ru + rd) / 4;
					for (int y_ind = 1; y_ind < nodes - 1; y_ind++)
					{
						for (int x_ind = 1; x_ind < nodes - 1; x_ind++)
						{
							vec3 point = gll::points[order - 1][x_ind] * vec_x
								+ gll::points[order - 1][y_ind] * vec_y
								+ mid;
							auto find = std::find_if(computational_mesh.nodes.begin(), computational_mesh.nodes.end(), [&](vec3 b) { return (point - b).norm() < 1e-16; });
							if (find == computational_mesh.nodes.end())
							{
								computational_mesh.elems.push_back(computational_mesh.nodes.size());
								computational_mesh.nodes.push_back(point);
							}
							else
							{
								computational_mesh.elems.push_back(find - computational_mesh.nodes.begin());
							}
						}
					}
				};

			for (auto& f : bound.faces)
			{
				add_face(f[0], f[1], f[2], f[3]);
			}

		}
		// volume
		{
			vec3 min   = original_mesh.nodes[original_mesh.elems[0 + original_offset]];
			vec3 x_max = original_mesh.nodes[original_mesh.elems[1 + original_offset]];
			vec3 y_max = original_mesh.nodes[original_mesh.elems[3 + original_offset]];
			vec3 z_max = original_mesh.nodes[original_mesh.elems[4 + original_offset]];

			vec3 vec_x = x_max - (x_max + min) / 2;
			vec3 vec_z = z_max - (z_max + min) / 2;
			vec3 vec_y = y_max - (y_max + min) / 2;

			vec3 mid_point = (
				original_mesh.nodes[original_mesh.elems[0 + original_offset]] +
				original_mesh.nodes[original_mesh.elems[1 + original_offset]] +
				original_mesh.nodes[original_mesh.elems[2 + original_offset]] +
				original_mesh.nodes[original_mesh.elems[3 + original_offset]] +
				original_mesh.nodes[original_mesh.elems[4 + original_offset]] +
				original_mesh.nodes[original_mesh.elems[5 + original_offset]] +
				original_mesh.nodes[original_mesh.elems[6 + original_offset]] +
				original_mesh.nodes[original_mesh.elems[7 + original_offset]]) / 8;

			for (int z_ind = 1; z_ind < nodes - 1; z_ind++)
			{
				for (int y_ind = 1; y_ind < nodes - 1; y_ind++)
				{
					for (int x_ind = 1; x_ind < nodes - 1; x_ind++)
					{
						vec3 point
							= gll::points[order - 1][x_ind] * vec_x
							+ gll::points[order - 1][y_ind] * vec_y
							+ gll::points[order - 1][z_ind] * vec_z
							+ mid_point;
						computational_mesh.elems.push_back(computational_mesh.nodes.size());
						computational_mesh.nodes.push_back(point);
					}
				}
			}
		}
	}

	fc::fc(std::filesystem::path file)
	{
		std::ifstream filestream(file, std::ios::in);
		if (!filestream)
		{
			throw std::runtime_error("cannot open fc file: " + file.string());
		}

		auto fc_file = nlohmann::json::parse(filestream);

		dim = 3;
		// dimensions
		if (fc_file["settings"]["dimensions"] == "2D") {
			dim = 2;
		}
		if (fc_file["settings"]["dimensions"] == "3D") {
			dim = 3;
		}

		// calc type
		if (fc_file["settings"]["type"] == "2D") {
			type = ploblem_type::eStatic;
		}
		if (fc_file["settings"]["type"] == "dynamic") {
			type = ploblem_type::eDynamic;
		}

		// materials
		unsigned int mat_id = 1;
		for (auto& mat : fc_file["materials"])
		{
			std::string e_b64 = mat["elasticity"][0]["constants"][0];
			std::string nu_b64 = mat["elasticity"][0]["constants"][1];

			double E, nu;
			base64::decode(e_b64.data(), e_b64.size(), reinterpret_cast<char*>(&E));
			base64::decode(nu_b64.data(), nu_b64.size(), reinterpret_cast<char*>(&nu));

			materials.emplace_back(mat_id, mat_id, E, nu);

			matid_threshold_map[mat["id"]] = mat_id;
			mat_id++;
		}

		// create block-threshold map
		for (auto& block : fc_file["blocks"])
		{
			block_threshold_map[block["id"]] = matid_threshold_map[block["material_id"]];
		}

		for (auto const& imap : block_threshold_map) {
			blocks.push_back(imap.first);
			thresholds.push_back(imap.second);
		}

		
		// mesh		

		size_t elems_count = fc_file["mesh"]["elems_count"];
		size_t nodes_count = fc_file["mesh"]["nodes_count"];
		original_mesh.elem_shifts.resize(elems_count + 1);
		base64::decode_vector(original_mesh.order, fc_file["mesh"]["elem_orders"]);
		base64::decode_vector(original_mesh.elem_type, fc_file["mesh"]["elem_types"]);
		base64::decode_vector(original_mesh.elemids, fc_file["mesh"]["elemids"]);
		base64::decode_vector(original_mesh.nids, fc_file["mesh"]["nids"]);
		base64::decode_vector(original_mesh.elems, fc_file["mesh"]["elems"]);		
		base64::decode_vector(original_mesh.nodes, fc_file["mesh"]["nodes"]);
		

		

		for (int node_id = 0; node_id < nodes_count; node_id++)
		{
			if (!original_mesh.map_node_numeration.insert_or_assign(original_mesh.nids[node_id], node_id).second)
			{
				throw std::runtime_error("Some nodes with the same ID in the mesh.\nToDo: Verify mesh nodes");
			}
		}
		for (int elem_id = 0; elem_id < elems_count; elem_id++) 
		{
			if (!original_mesh.map_element_numeration.insert_or_assign(original_mesh.elemids[elem_id], elem_id).second)
			{
				throw std::runtime_error("Some elements with the same ID in the mesh.\nToDo: Verify mesh elements");
			}
		}


		std::vector<int> elems_tmp = original_mesh.elems;
		int offset = 0;
		for (int elem_id = 0; elem_id < elems_count; elem_id++)
		{
			if (original_mesh.elem_type[elem_id] == '\x3')
			{
				// HEX (1st order)
				original_mesh.elem_shifts[elem_id] = offset;
				for (int i = offset; i < offset + 8; i++)
				{
					original_mesh.elems[i] = original_mesh.map_node_numeration[elems_tmp[i]];	
				}
				offset += 8;
				// cells with difficult cell types unsupported
				assert(!original_mesh.spectral_elems);
			}
			else if (original_mesh.elem_type[elem_id] == '\x4')
			{
				original_mesh.elem_shifts[elem_id] = offset;
				for (int i = offset; i < offset + 20; i++)
				{
					original_mesh.elems[i] = original_mesh.map_node_numeration[elems_tmp[i]];
				}
				offset += 20;
				original_mesh.spectral_elems = true;
			}
			else
			{
				// unsupported cell types
				assert(false);
			}
		}
		original_mesh.elem_shifts[elems_count] = offset;
		
		if(!original_mesh.spectral_elems)
		{
			computational_mesh = original_mesh;
		}
		else
		{
			//computational_mesh.elemids = original_mesh.elemids;
			computational_mesh.spectral_elems = true;
			computational_mesh.elem_type = original_mesh.elem_type;
			computational_mesh.order = original_mesh.order;
			computational_mesh.elem_shifts.resize(computational_mesh.elem_type.size() + 1);
			int offset_sem = 0;
			// remeshing
			for (int elem_id = 0; elem_id < elems_count; elem_id++)
			{
				computational_mesh.elem_shifts[elem_id] = offset_sem;
				assert(computational_mesh.elem_type[elem_id] == '\x4');
				add_spectral_elem(elem_id);
				int nodes = computational_mesh.order[elem_id] + 1;
				int nodes3 = nodes * nodes * nodes;
				offset_sem += nodes3;
			}
			computational_mesh.elem_shifts[elems_count] = offset_sem;
		}

		// loads
		for (auto& load : fc_file["loads"])
		{
			loads.push_back({});

			BC& bc = loads.back();
			base64::decode_vector(bc.apply_to, load["apply_to"]);
			bc.name = load["name"];
			int data_size;
			if (bc.name == "Pressure")
			{
				data_size = 1;
			}
			else if (bc.name == "Force")
			{
				bc.node_mapping(original_mesh, computational_mesh);
				data_size = 6;
			}
			else
			{
				throw std::runtime_error("not Force or Pressure load. Not supported yet");
			}
			for (int i = 0; i < data_size; i++)
			{
				std::string data_b64 = load["data"][i];
				bc.data.push_back(0.);
				bc.types.push_back({});
				if (type == ploblem_type::eDynamic)
				{
					std::string data = load["data"][i];
					std::string& name = bc.types[i].name;
					std::vector<double>& parameters = bc.types[i].param;
					auto parse = [&](const std::string& input) {

						// Извлекаем коэффициент
						std::stringstream preBracketStream(input);
						preBracketStream >> bc.data[i];

						size_t stBracketPos = input.find('(');
						size_t fiBracketPos = input.find_last_of(')');

						std::string dynamic = input.substr(stBracketPos + 1, fiBracketPos - stBracketPos - 1);
						
						stBracketPos = dynamic.find('(');
						fiBracketPos = dynamic.find_last_of(')');

						name = dynamic.substr(0, stBracketPos);

						std::stringstream paramsStream(dynamic.substr(stBracketPos + 1, fiBracketPos - stBracketPos - 1));
						double param;
						while (paramsStream >> param) {
							parameters.push_back(param);
							if (paramsStream.peek() == ',')
								paramsStream.ignore();
						}

					};
					parse(data);

				}
				else
				{
					base64::decode(data_b64.data(), data_b64.size(), reinterpret_cast<char*>(&bc.data[i]));
				}
				
				
			}

		}

		// restraints
		for (auto& restraint : fc_file["restraints"])
		{
			restraints.push_back({});

			BC& bc = restraints.back();
			base64::decode_vector(bc.apply_to, restraint["apply_to"]);
			int data_size = 6;
			for (int i = 0; i < data_size; i++)
			{
				std::string data_b64 = restraint["data"][i];
				bc.data.push_back(0.);
				base64::decode(data_b64.data(), data_b64.size(), reinterpret_cast<char*>(&bc.data[i]));

				bc.flag.push_back(int(restraint["flag"][i]));
			}
			bc.name = restraint["name"];

			bc.from_nodeset(original_mesh, computational_mesh);
		}

	}
}