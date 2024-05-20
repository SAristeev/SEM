#pragma once
#ifndef __FC_H__
#define __FC_H__

#include <num_types.h>

#include<filesystem>
#include<map>
#include<set>

namespace pre {

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
		material_t(unsigned char _id, int _threshold, double _E, double _nu) :id(_id), threshold(_threshold), E(_E), nu(_nu) {};
		unsigned char id;
		int threshold;

		double E;
		double nu;
	};

	struct UnstructedMesh {
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
	};

	struct BC {
		std::vector<int> apply_to;
		std::vector<bool> flag;
		std::vector<double> data;
		std::string name;
	};

	struct fc
	{
		int dim;

		std::map<int, unsigned char> matid_threshold_map;
		std::map<int, unsigned char> block_threshold_map;
		std::vector<material_t> materials;

		std::vector<int> blocks;
		std::vector<unsigned char> thresholds;

		UnstructedMesh mesh;

		std::vector<BC> restraints;
		std::vector<BC> loads;

		fc(std::filesystem::path path);
		
		void add_spectral_elem(int elem_id, int& offset, int& realoffset, const std::vector<int>& elems_tmp, const std::map<int, int>& map_node_numeration);
	};
};
#endif // __FC_H__
