#include "Common.h"
#include "Pipeline.h"

int main(int argc, char** argv)
{
    try
    {
        if(argc == 2 && std::string(argv[1]) == "--self-test")
        {
            runSelfTests();
            std::cout << "All A-BIRCH self-tests passed\n";
            return 0;
        }
        if(argc == 1) return runPipeline();
        if(argc == 10)
            return runPipeline(argv[1], argv[2], argv[3], argv[4],
                static_cast<size_t>(std::stoul(argv[5])),
                static_cast<size_t>(std::stoul(argv[6])),
                static_cast<size_t>(std::stoul(argv[7])),
                static_cast<size_t>(std::stoul(argv[8])), std::stod(argv[9]));
        throw std::invalid_argument(
            "Usage: birch [features truth results report branching min_cluster "
            "pilot_points pilot_max_k threshold_override]");
    }
    catch(const std::exception& error)
    {
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }
}
