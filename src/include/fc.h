#pragma once
#ifndef __FC_H__
#define __FC_H__

#include<Eigen/Core>
#include<filesystem>
#include<map>

namespace pre {

	using vec3 = Eigen::Vector3d;

	struct material_t
	{
		material_t(unsigned char _id, int _threshold, double _E, double _nu) :id(_id), threshold(_threshold), E(_E), nu(_nu) {};
		unsigned char id;
		int threshold;

		double E;
		double nu;
	};

	struct UnstructedMesh {
		std::vector<vec3> nodes;
		std::vector<int>          nids;
		std::vector<int>          elems;
		std::vector<int>          elemids;
		std::vector<int>          elem_shifts;
		std::vector<int>          order;
		std::vector<uint8_t>      elem_type;

		std::unordered_map<int, std::tuple<int, int, int>> hex_loc_indxs;
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

		void set_order();
	};
};
#endif // __FC_H__
