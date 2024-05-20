#include <fc.h>
#include <solver.h>


#include <unordered_map>
#include <string>

using namespace pre;
using namespace solver;

int main(int argc, char* argv[]) {    
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

    fc fcase(FC_FileName);
    solve(fcase, VTU_FileName);
    
    return 0;
}


