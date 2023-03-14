#include "qx/Simulator.hpp"

#include "qx/DefaultOperations.hpp"
#include "qx/LibqasmInterface.hpp"
#include "qx/Version.hpp"

namespace {

struct CLIOptions {
    std::uint64_t iterations = 0;
    std::string filePath;
};

std::optional<CLIOptions> parseCLIOptions(int argc, char** argv) {
    if (argc == 2) {
        return {{.iterations = 0, .filePath = std::string(argv[1])}};
    }

    if (argc == 4 && std::string(argv[1]) == "-c") {
        auto iterations = atoi(argv[2]);

        if (iterations < 0) {
            return std::nullopt;
        }

        return {{.iterations = static_cast<std::uint64_t>(iterations), .filePath = std::string(argv[3])}};
    }

    return std::nullopt;
}

void printUsage(char** argv) {
    std::cerr << "Usage: \n   " << argv[0] << " [-c iterations] <file>\n"
                  << "\titerations: the number of shots, >= 0\n"
                  << "\nfile: the cQasm file to execute\n"
                  << std::endl;
}

}

void print_banner() {
    // clang-format off
    std::cout << "\n=================================================================================================== \n"; 
    std::cout << "        _______                                                                                       \n";
    std::cout << "       /  ___   \\   _  __      ____   ____   __  ___  __  __   __    ___  ______  ____    ___        \n";
    std::cout << "      /  /   /  |  | |/ /     / __/  /  _/  /  |/  / / / / /  / /   / _ |/_  __/ / __ \\  / _ \\      \n";
    std::cout << "     /  /___/  /   >   <     _\\ \\   _/ /   / /|_/ / / /_/ /  / /__ / __ | / /   / /_/ / / , _/      \n";
    std::cout << "     \\______/\\__\\ /_/|_|    /___/  /___/  /_/  /_/  \\____/  /____//_/ |_|/_/    \\____/ /_/|_|    \n";
    std::cout << "                                                                                                      \n";
    std::cout << "       Version " << QX_VERSION << " - QuTech - " << QX_RELEASE_YEAR << " - report bugs and suggestions to: p.lehenaff@tudelft.nl\n";
    std::cout << "  =================================================================================================== \n";
    std::cout << "\n" << std::endl;
    // clang-format on
}

int main(int argc, char **argv) {
    print_banner();

    auto cliOptions = parseCLIOptions(argc, argv);

    if (!cliOptions) {
        printUsage(argv);
        return -1;
    }

    std::cout << "Will produce samples for " << cliOptions->iterations << " iteration"
              << (cliOptions->iterations > 1 ? "s" : "") << " from cQasm file '" << cliOptions->filePath
              << "'..." << std::endl;

    auto output = qx::executeFile(
        cliOptions->filePath, qx::default_operations::defaultOperations, cliOptions->iterations);

    if (auto *simulationError = std::get_if<qx::SimulationError>(&output)) {
        std::cerr << "Simulation failed:\n"
                  << simulationError->message << std::endl;
        return 1;
    }

    std::cout << std::get<qx::SimulationResult>(output) << std::endl;
    return 0;
}