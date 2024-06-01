#pragma once
#ifndef __VTK_H__
#define __VTK_H__

#include "fc.h"

#include <vector>
#include <filesystem>
namespace post
{
	void export2vtk(const pre::UnstructedMesh& fcase, const std::vector<double>& results, const std::filesystem::path& filename);
	void collect_steps(const std::filesystem::path& dir, const std::string& filename, const std::vector<double>& time_steps);
}
#endif __VTK_H__
