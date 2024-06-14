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
	void export2vtk(const pre::UnstructedMesh& mesh, const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& a, const std::vector<std::vector<double>>& eps, const std::vector<std::vector<double>>& sigma, const std::vector<double>& F, const std::filesystem::path& filename)
	{
		vtkNew<vtkUnstructuredGrid> unstructuredGrid;

		// for debug
		vtkIndent ind;

		// points
		vtkNew<vtkPoints> points;
		vtkNew<vtkDoubleArray> u_array;
		vtkNew<vtkDoubleArray> v_array;
		vtkNew<vtkDoubleArray> a_array;

		vtkNew<vtkDoubleArray> eps_array;
		vtkNew<vtkDoubleArray> sigma_array;

		vtkNew<vtkDoubleArray> F_array;

		{
			u_array->SetName("Displacement");
			u_array->SetNumberOfComponents(3);
			u_array->SetNumberOfTuples(mesh.nodes.size());

			v_array->SetName("Velocity");
			v_array->SetNumberOfComponents(3);
			v_array->SetNumberOfTuples(mesh.nodes.size());

			a_array->SetName("Acceleration");
			a_array->SetNumberOfComponents(3);
			a_array->SetNumberOfTuples(mesh.nodes.size());
		}
		{
			eps_array->SetName("Strain");
			eps_array->SetNumberOfComponents(6);
			eps_array->SetNumberOfTuples(mesh.nodes.size());
		
			sigma_array->SetName("Stress");
			sigma_array->SetNumberOfComponents(6);
			sigma_array->SetNumberOfTuples(mesh.nodes.size());
		}


		F_array->SetName("F");
		F_array->SetNumberOfComponents(3);
		F_array->SetNumberOfTuples(mesh.nodes.size());

		for (int node = 0; node < mesh.nodes.size(); node++)
		{
			points->InsertNextPoint(mesh.nodes[node][0], mesh.nodes[node][1], mesh.nodes[node][2]);
			double tuple_u[3] = { u[3 * node], u[3 * node + 1], u[3 * node + 2] };
			double tuple_v[3] = { v[3 * node], v[3 * node + 1], v[3 * node + 2] };
			double tuple_a[3] = { a[3 * node], a[3 * node + 1], a[3 * node + 2] };
			double tuple_F[3] = { F[3 * node], F[3 * node + 1], F[3 * node + 2] };
			u_array->SetTuple(node, tuple_u);
			v_array->SetTuple(node, tuple_v);
			a_array->SetTuple(node, tuple_a);
			F_array->SetTuple(node, tuple_F);


			double tuple_eps[6] = { eps[0][node], eps[1][node], eps[2][node], eps[3][node], eps[4][node], eps[5][node] };
			double tuple_sigma[6] = { sigma[0][node], sigma[1][node], sigma[2][node], sigma[3][node], sigma[4][node], sigma[5][node] };
			eps_array->SetTuple(node, tuple_eps);
			sigma_array->SetTuple(node, tuple_sigma);
		}
	

		unstructuredGrid->SetPoints(points);
		unstructuredGrid->GetPointData()->AddArray(u_array);
		unstructuredGrid->GetPointData()->AddArray(v_array);
		unstructuredGrid->GetPointData()->AddArray(a_array);
		unstructuredGrid->GetPointData()->AddArray(F_array);

		unstructuredGrid->GetPointData()->AddArray(eps_array);
		unstructuredGrid->GetPointData()->AddArray(sigma_array);

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

	void collect_steps(const std::filesystem::path& dir, const std::string& filename, const std::vector<int>& int_steps, const std::vector<double>& time_steps)
	{
		std::ofstream pvd_file;
		std::string full_name = filename + ".pvd";
		pvd_file.open(dir / full_name);
		pvd_file << "<?xml version=\"1.0\"?>\n";
		pvd_file << " <VTKFile type=\"Collection\" version=\"0.1\">\n";
		pvd_file << "  <Collection>\n";
		pvd_file << "  <DataSet timestep= \"0\" step=\"0\" file=\""	<< std::string("results/" + filename + "_0.vtu\"/>\n");
		for (int i = 1; i < time_steps.size(); i++) 
		{
			pvd_file << "  <DataSet timestep=\""
				<< time_steps[i] << "\" step=\""
				<< i
				<< "\" file=\""
				<< std::string("results/" + filename + "_" + std::to_string(int_steps[i]) + ".vtu\"/>\n");

		};
		pvd_file << "  </Collection>\n";
		pvd_file << "</VTKFile>\n";
		pvd_file.close();
	}
}