#pragma once
#ifndef __FC_H__
#define __FC_H__

#include <num_types.h>

#include <functional>
#include <filesystem>
#include <map>
#include <set>

namespace pre {
	enum ploblem_type {
		eStatic,
		eDynamic
	};
	using namespace types;

	struct parametic_points
	{
		std::vector<vec3> v;
		vec3 mid = { 0.5,0.5,0.5 };
		parametic_points()
		{
			v.resize(8);
			v[0] = { 0,0,0 };
			v[1] = { 1,0,0 };
			v[2] = { 1,1,0 };
			v[3] = { 0,1,0 };
			v[4] = { 0,0,1 };
			v[5] = { 1,0,1 };
			v[6] = { 1,1,1 };
			v[7] = { 0,1,1 };
		}
	};

	struct material_t
	{
		material_t(unsigned char _id, int _threshold, double _E, double _nu, double _density) :id(_id), threshold(_threshold), E(_E), nu(_nu), density(_density) {};
		unsigned char id;
		int threshold;

		double E;
		double nu;
		double density;
	};

	//TODO: Add reverse Cuthill–McKee
	struct UnstructedMesh {
		bool spectral_elems;
		// array of all points
		std::vector<vec3>         nodes;
		std::vector<int>          nids;

		// connectivity analog
		// contains node's id by elems
		std::vector<int>          elems;
		// shifts of elems
		std::vector<int>          elem_shifts;

		// elems parameters 
		std::vector<int>          elemids;
		std::vector<uint8_t>      elem_type;
		std::vector<int>          order;	

		std::map<int, int>        map_node_numeration;
		std::map<int, int>        map_element_numeration;
	};

	// TODO: interpolator
	struct BC {
		struct dynamic_type
		{
			std::string name;
			std::vector<double> param;
		}; 

		std::vector<dynamic_type> types;
		std::vector<int> apply_to;
		std::vector<bool> flag;
		std::vector<double> data;

		std::string name;

		std::vector<std::pair<int, int>> elem2face;
		void from_nodeset(const UnstructedMesh& original_mesh, const UnstructedMesh& computational_mesh);
		void node_mapping(const UnstructedMesh& original_mesh, const UnstructedMesh& computational_mesh);

		
		
	};

	struct fc
	{
		ploblem_type type; // 0 static, 1 dymamic 
		int dim;
		std::vector<std::map<int, int>> elem_face2_nodes;
		std::map<int, unsigned char> matid_threshold_map;
		std::map<int, unsigned char> block_threshold_map;
		std::vector<material_t> materials;

		std::vector<int> blocks;
		std::vector<unsigned char> thresholds;

		UnstructedMesh original_mesh;
		UnstructedMesh computational_mesh;

		std::vector<BC> restraints;
		std::vector<BC> loads;

		fc(std::filesystem::path path);
	private:
		void add_spectral_elem(int elem_id);
	};
};
#endif // __FC_H__
