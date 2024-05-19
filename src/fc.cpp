#include <fc.h>

#include <fstream>
#include <algorithm> 

#include <base64.h>
#include <json/json.hpp>

#include <vtkCell.h>

namespace pre {
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
		base64::decode_vector(mesh.nodes, fc_file["mesh"]["nodes"]);

		std::map<int, int>        map_node_numeration;
		std::map<int, int>        map_element_numeration;

		for (int node_id = 0; node_id < nodes_count; node_id++)
		{
			if (!map_node_numeration.insert_or_assign(mesh.nids[node_id], node_id).second)
			{
				throw std::runtime_error("Some nodes with the same ID in the mesh.\nToDo: Verify mesh nodes");
			}
		}

		// Only hex mesh
		int offset = 0;
		for (int elem_ID = 0; elem_ID < elems_count; elem_ID++)
		{
			if (!map_element_numeration.insert_or_assign(mesh.elemids[elem_ID], elem_ID).second)
			{
				throw std::runtime_error("Some elements with the same ID in the mesh.\nToDo: Verify mesh elements");
			}
			if (mesh.elem_type[elem_ID] == '\x3') 
			{
				// HEX (1st order)
				mesh.elem_shifts[elem_ID] = offset;
				
				for (int i = offset; i < offset + 8; i++)
				{
					mesh.elems[i] = map_node_numeration[elems_tmp[i]];
				}
				offset += 8;
			}
			else if (mesh.elem_type[elem_ID] == '\x4')
			{
				// HEX (2nd and higher order)
				mesh.elem_shifts[elem_ID] = offset;
				for (int i = offset; i < offset + 20; i++)
				{
					mesh.elems[i] = map_node_numeration[elems_tmp[i]];
				}
				offset += 20;
			}

		}
		mesh.elem_shifts[elems_count] = offset;

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