#include "fc.h"
#include "solver.h"

#include <filesystem>
#include <unordered_map>
#include <string>

using namespace pre;
using namespace solver;


//TODO: make it parallel (shared memory)
// maybe CUDA-powered algebra
// or MPI version

int main(int argc, char* argv[]) {    
    std::unordered_map<std::string, std::string> parsed_params;//in the pair {key,param} param may be empty

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
    
    std::filesystem::path WorkDirectory;
    {
        const auto it = parsed_params.find("-d");
        if (it != parsed_params.end() && !it->second.empty()) {
            WorkDirectory = it->second;
            parsed_params.erase(it);
        }
        else {
            printf("ERROR: Path to fc is not provided!\n");
        }
    }

    std::string filename;
    {
        const auto it = parsed_params.find("-f");
        if (it != parsed_params.end() && !it->second.empty()) {
            filename = it->second;
            parsed_params.erase(it);
        }
        else {
            printf("ERROR: fc filename is not provided!\n");
        }
    }

    fc fcase(WorkDirectory / std::string(filename + +".fc"));
    start_problem(fcase, WorkDirectory, filename );
    
    return 0;
}


