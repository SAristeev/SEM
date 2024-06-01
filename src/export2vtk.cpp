#include "export2vtk.h"

#include <fstream>
#include <vector>

#include <vtkPointData.h>
#include <vtkDoubleArray.h>
#include <vtkCellArray.h>
#include <vtkNew.h>
#include <vtkHexahedron.h>
#include <vtkQuadraticHexahedron.h>
#include <vtkLagrangeHexahedron.h>

#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>

namespace post {
	void export2vtk(const pre::UnstructedMesh& mesh,const std::vector<double>& results, const std::filesystem::path& filename)
	{		
		vtkNew<vtkUnstructuredGrid> unstructuredGrid;
		
		// for debug
		vtkIndent ind;

		// points
		vtkNew<vtkPoints> points;
		vtkNew<vtkDoubleArray> dataArray;
		dataArray->SetName("Displacement");
		dataArray->SetNumberOfComponents(3);
		dataArray->SetNumberOfTuples(mesh.nodes.size());

		for (int node = 0; node < mesh.nodes.size(); node++)
		{
			points->InsertNextPoint(mesh.nodes[node][0], mesh.nodes[node][1], mesh.nodes[node][2]);
			double tuple[3] = { results[3 * node], results[3 * node + 1], results[3 * node + 2] };
			dataArray->SetTuple(node, tuple);
		}
	

		unstructuredGrid->SetPoints(points);
		unstructuredGrid->GetPointData()->AddArray(dataArray);

		// elems
		 
		// all cells MUST have equal type
		VTKCellType type;
		vtkNew<vtkCellArray> cellArray;
		int offset = 0;
		for (int elem_id = 0; elem_id < mesh.elem_type.size(); elem_id++)
		{
			if (mesh.elem_type[elem_id] == '\x3')
			{
				type = VTK_HEXAHEDRON;
				vtkNew<vtkHexahedron> hex;
				for (int i = 0; i < 8; i++)
				{
					hex->GetPointIds()->SetId(i, mesh.elems[i + offset]);
				}
				offset += 8;
				cellArray->InsertNextCell(hex);
			}
			if (mesh.elem_type[elem_id] == '\x4')
			{
				type = VTK_LAGRANGE_HEXAHEDRON;
				vtkNew<vtkLagrangeHexahedron> hex;
				int nodes = mesh.order[elem_id] + 1;
				int nodes3 = nodes * nodes * nodes;
					
				hex->GetPointIds()->SetNumberOfIds(nodes3);
				for (int i = 0; i < nodes3; i++)
				{
					hex->GetPointIds()->SetId(i, mesh.elems[i + offset]);
				}

				cellArray->InsertNextCell(hex);
				offset += nodes3;
			}
		}
		
		unstructuredGrid->SetCells(type, cellArray);

		// Write file.
		vtkNew<vtkXMLUnstructuredGridWriter> writer;
		writer->SetFileName(filename.string().c_str());
		writer->SetInputData(unstructuredGrid);
		writer->Write();

	}
	void collect_steps(const std::filesystem::path& dir, const std::string& filename, const std::vector<double>& time_steps)
	{
		std::ofstream pvd_file;
		std::string full_name = filename + ".pvd";
		pvd_file.open(dir / full_name);
		pvd_file << "<?xml version=\"1.0\"?>\n";
		pvd_file << " <VTKFile type=\"Collection\" version=\"0.1\">\n";
		pvd_file << "  <Collection>\n";
		pvd_file << "  <DataSet timestep= \"0\" step=\"0\" file=\""	<< std::string(filename + "_0.vtu\"/>\n");
		for (int i = 1; i < time_steps.size(); i++) 
		{
			pvd_file << "  <DataSet timestep=\""
				<< time_steps[i] << "\" step=\""
				<< i
				<< "\" file=\""
				<< std::string(filename + "_" + std::to_string(i) + ".vtu\"/>\n");

		};
		pvd_file << "  </Collection>\n";
		pvd_file << "</VTKFile>\n";
		pvd_file.close();
	}
}