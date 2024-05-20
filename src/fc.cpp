#include <fc.h>
#include <gll.h>
#include <parametric_hex.h>

#include <fstream>
#include <algorithm>	

#include <base64.h>
#include <json/json.hpp>

namespace pre {
	static boundary_indexes bound;
	void fc::add_spectral_elem(int elem_id, int offset_old, int offset_real, const std::vector<int>& elems_tmp, const std::vector<vec3>& nodes_tmp, const std::map<int, int>& map_node_numeration)
	{
		int order = mesh.order[elem_id];
		int nodes = order + 1;
		int nodes3 = nodes * nodes * nodes;
		mesh.elems.resize(mesh.elems.size() + nodes3 - 20);

		// vertices
		{
			for (int i = offset_old; i < offset_old + 8; i++)
			{
				int node = map_node_numeration.at(elems_tmp[i]);
				vec3 point = nodes_tmp[node];
				auto find = std::find_if(mesh.nodes.begin(), mesh.nodes.end(), [&](vec3 b) { return (point - b).norm() < 1e-16; });
				if (find == mesh.nodes.end())
				{
					mesh.elems[offset_real] = mesh.nodes.size();
					mesh.nodes.push_back(point);
					offset_real++;
					mesh.real_nodes++;
				}
				else
				{
					mesh.elems[offset_real] = find - mesh.nodes.begin();
					offset_real++;
				}
			}
		}
		// edges
		{
			auto add_edge = [&](int x_min_idx, int x_max_idx)
				{
					vec3 x_min = nodes_tmp[map_node_numeration.at(elems_tmp[x_min_idx + offset_old])];
					vec3 x_max = nodes_tmp[map_node_numeration.at(elems_tmp[x_max_idx + offset_old])];
					vec3 x_mid = (x_max + x_min) / 2;
					vec3 x_vec = x_max - x_mid;

					for (int i = 1; i < nodes - 1; i++)
					{
						vec3 point = gll::points[order - 1][i] * x_vec + x_mid;
						auto find = std::find_if(mesh.nodes.begin(), mesh.nodes.end(), [&](vec3 b) { return (point - b).norm() < 1e-16; });
						if (find == mesh.nodes.end())
						{
							mesh.elems[offset_real] = mesh.nodes.size();
							mesh.nodes.push_back(point);
							offset_real++;
							mesh.real_nodes++;
						}
						else
						{
							mesh.elems[offset_real] = find - mesh.nodes.begin();
							offset_real++;
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
					vec3 ld = nodes_tmp[map_node_numeration.at(elems_tmp[ld_indx + offset_old])];
					vec3 rd = nodes_tmp[map_node_numeration.at(elems_tmp[rd_indx + offset_old])];
					vec3 lu = nodes_tmp[map_node_numeration.at(elems_tmp[lu_indx + offset_old])];
					vec3 ru = nodes_tmp[map_node_numeration.at(elems_tmp[ru_indx + offset_old])];

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
							auto find = std::find_if(mesh.nodes.begin(), mesh.nodes.end(), [&](vec3 b) { return (point - b).norm() < 1e-16; });
							if (find == mesh.nodes.end())
							{
								mesh.elems[offset_real] = mesh.nodes.size();
								mesh.nodes.push_back(point);
								mesh.real_nodes++;
								offset_real++;
							}
							else
							{
								mesh.elems[offset_real++] = find - mesh.nodes.begin();
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
			vec3 min   = nodes_tmp[map_node_numeration.at(elems_tmp[0 + offset_old])];
			vec3 x_max = nodes_tmp[map_node_numeration.at(elems_tmp[1 + offset_old])];
			vec3 y_max = nodes_tmp[map_node_numeration.at(elems_tmp[3 + offset_old])];
			vec3 z_max = nodes_tmp[map_node_numeration.at(elems_tmp[4 + offset_old])];

			vec3 vec_x = x_max - (x_max + min) / 2;
			vec3 vec_z = z_max - (z_max + min) / 2;
			vec3 vec_y = y_max - (y_max + min) / 2;

			vec3 mid_point = (
				nodes_tmp[map_node_numeration.at(elems_tmp[0 + offset_old])] +
				nodes_tmp[map_node_numeration.at(elems_tmp[1 + offset_old])] +
				nodes_tmp[map_node_numeration.at(elems_tmp[2 + offset_old])] +
				nodes_tmp[map_node_numeration.at(elems_tmp[3 + offset_old])] +
				nodes_tmp[map_node_numeration.at(elems_tmp[4 + offset_old])] +
				nodes_tmp[map_node_numeration.at(elems_tmp[5 + offset_old])] +
				nodes_tmp[map_node_numeration.at(elems_tmp[6 + offset_old])] +
				nodes_tmp[map_node_numeration.at(elems_tmp[7 + offset_old])]) / 8;

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

						mesh.elems[offset_real++] = mesh.nodes.size();
						mesh.nodes.push_back(point);
						mesh.real_nodes++;
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
		mesh.elem_shifts.resize(elems_count + 1);

		base64::decode_vector(mesh.order, fc_file["mesh"]["elem_orders"]);
		base64::decode_vector(mesh.elem_type, fc_file["mesh"]["elem_types"]);
		base64::decode_vector(mesh.elemids, fc_file["mesh"]["elemids"]);
		base64::decode_vector(mesh.nids, fc_file["mesh"]["nids"]);

		std::vector<int> elems_tmp;
		base64::decode_vector(elems_tmp, fc_file["mesh"]["elems"]);
		mesh.elems.resize(elems_tmp.size());
		
		std::vector<vec3> nodes_tmp;
		base64::decode_vector(nodes_tmp, fc_file["mesh"]["nodes"]);
		

		std::map<int, int>        map_node_numeration;
		std::map<int, int>        map_element_numeration;

		for (int node_id = 0; node_id < nodes_count; node_id++)
		{
			if (!map_node_numeration.insert_or_assign(mesh.nids[node_id], node_id).second)
			{
				throw std::runtime_error("Some nodes with the same ID in the mesh.\nToDo: Verify mesh nodes");
			}
		}
		for (int elem_id = 0; elem_id < elems_count; elem_id++) 
		{
			if (!map_element_numeration.insert_or_assign(mesh.elemids[elem_id], elem_id).second)
			{
				throw std::runtime_error("Some elements with the same ID in the mesh.\nToDo: Verify mesh elements");
			}
		}
		bool spectral_elems = false;
		int offset_old = 0;
		int offset_real = 0;
		for (int elem_id = 0; elem_id < elems_count; elem_id++)
		{
			if (mesh.elem_type[elem_id] == '\x3')
			{
				// HEX (1st order)
				mesh.elem_shifts[elem_id] = offset_real;
				
				for (int i = offset_old; i < offset_old + 8; i++)
				{
					int node = map_node_numeration[elems_tmp[i]];
					mesh.elems[i] = mesh.nodes.size();
					mesh.nodes.push_back(nodes_tmp[node]);
				}
				offset_old += 8;
				offset_real += 8;
			}
			else if (mesh.elem_type[elem_id] == '\x4')
			{
				// HEX (2nd and higher order)
				mesh.elem_shifts[elem_id] = offset_real;
				
				if (mesh.order[elem_id] == 2) {
					for (int i = offset_old; i < offset_old + 20; i++)
					{
						int node = map_node_numeration[elems_tmp[i]];
						mesh.elems[i] = mesh.nodes.size();
						mesh.nodes.push_back(nodes_tmp[node]);
					}	
				}
				else
				{
					add_spectral_elem(elem_id, offset_old, offset_real, elems_tmp, nodes_tmp, map_node_numeration);
					spectral_elems = true;
				}
				offset_old += 20;
				offset_real += (mesh.order[elem_id] + 1) * (mesh.order[elem_id] + 1) * (mesh.order[elem_id] + 1);
			}
		}
		mesh.elem_shifts[elems_count] = offset_real;
		
		if (!spectral_elems) 
		{
			mesh.real_nodes = mesh.nodes.size();
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
				base64::decode(data_b64.data(), data_b64.size(), reinterpret_cast<char*>(&bc.data[i]));
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
		}
	}
}