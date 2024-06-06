#pragma once
#ifndef __VTK_H__
#define __VTK_H__

#include "fc.h"

#include <vector>
#include <filesystem>
namespace post
{
	void export2vtk(const pre::UnstructedMesh& mesh, const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& a, const std::vector<double>& F, const std::filesystem::path& filename);
	void collect_steps(const std::filesystem::path& dir, const std::string& filename, const std::vector<int>& int_steps, const std::vector<double>& time_steps);
}
#endif __VTK_H__
