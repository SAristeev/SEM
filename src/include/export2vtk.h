#pragma once
#ifndef __VTK_H__
#define __VTK_H__

#include <fc.h>
#include <vector>
#include <filesystem>
namespace post
{
	void export2vtk(const pre::UnstructedMesh& fcase, std::vector<double>& results, std::filesystem::path filename);
	void collect_steps(std::filesystem::path dir, std::string filename, std::vector<double> steps);
}
#endif __VTK_H__
