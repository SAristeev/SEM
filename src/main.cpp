
#include <tbb/tick_count.h>
#include <tbb/parallel_for.h>



#include <vtkActor.h>
#include <vtkDataSetMapper.h>
#include <vtkNamedColors.h>
#include <vtkNew.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkXMLUnstructuredGridReader.h>


#include <vtkActor.h>
#include <vtkCellArray.h>
#include <vtkDataSetMapper.h>
#include <vtkNamedColors.h>
#include <vtkNew.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkTetra.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include <vtkDoubleArray.h>
#include <vtkHigherOrderHexahedron.h>
#include <vtkHigherOrderQuadrilateral.h>
#include <vtkLagrangeQuadrilateral.h>
#include <vtkLagrangeHexahedron.h>

int main(int argc, char* argv[]) {
    tbb::tick_count start = tbb::tick_count::now();
    std::unordered_map<std::string, std::string>  parsed_params;//in the pair {key,param} param may be empty

    for (int pos = 1; pos < argc; ++pos) {
        if (argv[pos][0] == '-') {//key is found
            std::string key(argv[pos]);

            if (pos + 1 < argc && argv[pos + 1][0] != '-') {//value is found
                std::string value(argv[pos + 1]);

                parsed_params.insert(std::make_pair(key, value));
                ++pos;
                continue;
            }
            parsed_params.insert(std::make_pair(key, std::string()));
        }
    }
    
    std::string FC_FileName;
    {
        const auto it = parsed_params.find("--fc");
        if (it != parsed_params.end() && !it->second.empty()) {
            FC_FileName = it->second;
            parsed_params.erase(it);
        }
        else {
            printf("ERROR: Path to fc is not provided!\n");
        }
    }

    std::string VTU_FileName;
    {
        const auto it = parsed_params.find("--vtu");
        if (it != parsed_params.end() && !it->second.empty()) {
            VTU_FileName = it->second;
            parsed_params.erase(it);
        }
        else {
            printf("ERROR: Path to vtu is not provided!\n");
        }
    }


    return 0;
}


