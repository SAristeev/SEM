#include <export2vtk.h>

#include <iostream>
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
	void export2vtk(const pre::UnstructedMesh& mesh, std::vector<double>& results, std::filesystem::path filename)
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
				else 
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
		}
		
		unstructuredGrid->SetCells(type, cellArray);

		// Write file.
		vtkNew<vtkXMLUnstructuredGridWriter> writer;
		writer->SetFileName(filename.string().c_str());
		writer->SetInputData(unstructuredGrid);
		writer->Write();

	}
}