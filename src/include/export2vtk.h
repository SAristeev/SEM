#pragma once
#ifndef __VTK_H__
#define __VTK_H__

#include <fc.h>

#include <filesystem>
namespace post
{
	void export2vtk(const pre::UnstructedMesh& fcase, std::filesystem::path filename);
}
#endif __VTK_H__
