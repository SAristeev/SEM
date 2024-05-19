#include <export2vtk.h>

#include <iostream>

#include <vtkCellArray.h>
#include <vtkNew.h>
#include <vtkHexahedron.h>
#include <vtkQuadraticHexahedron.h>

#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>

namespace post {
	void export2vtk(const pre::UnstructedMesh& mesh, std::filesystem::path filename)
	{		
		vtkNew<vtkUnstructuredGrid> unstructuredGrid;
		
		// for debug
		vtkIndent ind;

		// points
		vtkNew<vtkPoints> points;
		for (int node = 0; node < mesh.nodes.size(); node++)
		{
			points->InsertNextPoint(mesh.nodes[node][0], mesh.nodes[node][1], mesh.nodes[node][2]);
		}
		unstructuredGrid->SetPoints(points);
				
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
				if(mesh.order[elem_id] == 2){
					type = VTK_QUADRATIC_HEXAHEDRON;
					vtkNew<vtkQuadraticHexahedron> hex;
					for (int i = 0; i < 20; i++)
					{
						hex->GetPointIds()->SetId(i, mesh.elems[i + offset]);
					}
					cellArray->InsertNextCell(hex);
					offset += 20;
				}
			}
		}
		
		unstructuredGrid->SetCells(type, cellArray);

		// Write file.
		vtkNew<vtkXMLUnstructuredGridWriter> writer;
		writer->SetFileName(filename.string().c_str());
		writer->SetInputData(unstructuredGrid);
		writer->Write();

	}
}