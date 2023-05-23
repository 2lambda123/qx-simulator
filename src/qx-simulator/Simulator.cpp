#include "qx/Simulator.hpp"

#include "qx/LibqasmInterface.hpp"
#include "qx/Version.hpp"
#include "qx/DefaultOperations.hpp"

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
    std::string filePath = "";
    size_t iterations = 1;
    print_banner();

    // I would put this in a different function
    // Maybe returning a std::expected<CommandLineParams>?
    // Where CommandLineParams contains a file and an optional iterations?
    int argIndex = 1;
    bool argParsingFailed = false;
    while (argIndex < argc) {
        auto currentArg = argv[argIndex];

        if (std::string(currentArg) == "-c") {
            if (argIndex + 1 >= argc) {
                argParsingFailed = true;
            } else {
                iterations = atoi(argv[++argIndex]);
            }
        } else {
            if (argIndex + 1 < argc) {
                argParsingFailed = true;
            } else {
                filePath = std::string(currentArg);
            }
        }

        if (argParsingFailed) {
            break;
        }

        ++argIndex;
    }

    if (filePath == "" || argParsingFailed) {
        // I would put the usage in a different function
        // And add some more information
        // E.g.
        // - a generic command line, without argv[0] (this can be an absolute path, can't it?)
        //   qx-simulator [-c <ITERATIONS>] <INPUT FILE>
        // - a description of each argument
        //   ITERATIONS: optional. Description of what these iterations mean
        //   INPUT FILE: CQASM assembler file.
        // - and some examples
        //   qx-simulator bell.qc
        //   qx-simulator -c 2500 bell.qc
        std::cerr << "Usage:\n   "
                  << argv[0] << " [-c iterations] file.qc" << std::endl;
        return -1;
    }

    // I would put this in a different function as well
    std::cout << "Will produce samples for "
              << iterations << " iteration" << (iterations > 1 ? "s" : "")
              << " from cQasm file '" << filePath << "'..." << std::endl;

    auto output = qx::executeFile(filePath, qx::default_operations::defaultOperations, iterations);

    if (auto* simulationError = std::get_if<qx::SimulationError>(&output)) {
        std::cerr << "Simulation failed:\n" << simulationError->message << std::endl;
        return 1;
    }

    std::cout << std::get<qx::SimulationResult>(output) << std::endl;
    return 0;
}

/*
struct CommandLineParams {
    std::string filePath;
    std::optional<size_t> iterations;
}

void print_usage() {}

int executeFile(std::string filePath, size_t iterations) {}

int main(int argc, char **argv) {
    print_banner();

    if (const auto result = parseCommandLine(argc, argv); result.has_value()) {
        auto clParams = *result;
        return executeFile(clParams.filePath, clParams.iterations);
    } else {
        print_usage();
        return -1;
    }
    return 0;
}
*/