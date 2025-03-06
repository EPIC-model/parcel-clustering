#include <iostream>
#include <string>
#include <tuple>
#include <list>
#include <cstdlib>
#include <filesystem>

namespace fs = std::filesystem;

using strint_t = std::tuple<std::string, std::string, int>;
using strlist_t = std::tuple<std::string, std::string, std::list<int>>;
using strbool_t = std::tuple<std::string, std::string, bool>;
using strstr_t = std::tuple<std::string, std::string, std::string>;

struct args_t {
    strint_t nSamples       = {"--nsamples", "number of random samples", 10};
    strlist_t nRanks        = {"--nranks", "number of MPI ranks for the parallel runs", {}};
    strint_t nParcelPerCell = {"--nppc", "number of parels (on average) per grid cell", 40};
    strint_t minVratio      = {"--min-vratio", "minimum volume ratio", 40.0};
    strint_t nx             = {"--nx", "number of grid cells in x", 32};
    strint_t ny             = {"--ny", "number of grid cells in y", 32};
    strint_t nz             = {"--nz", "number of grid cells in z", 32};
    strint_t seed           = {"--seed", "rng seed", 42};
    strint_t nTasksPerNode  = {"--ntasks-per-node", "number of tasks per node", -1};
    strstr_t commType       = {"--comm-type", "communication layer for the parallel version",  "p2p"};
    strbool_t shuffle       = {"--shuffle", "shuffle initial parcel locally", false};
    strbool_t subcomm       = {"--subcomm", "use subcommunicator in resolve tree step", false};
    strbool_t verbose       = {"--verbose", "print intermediate output", false};
    strstr_t cmd            = {"--cmd", "run jobs with 'mpirun' or 'srun'", "srun"};
};

void check_single_argument(char* argv[], int i) {
    if (argv[i+1] == nullptr || std::string(argv[i+1]).contains("--")) {
        throw std::invalid_argument("No value provided for " + std::string(argv[i]));
    }
}

int check_multiple_argument(char* argv[], int i) {
    int result = 0;

    int j = i+1;
    while (argv[j] != nullptr && !std::string(argv[j]).contains("--")) {
        j++;
        result++;
    }

    if (result == 0)  {
        throw std::invalid_argument("No value provided for " + std::string(argv[i]));
    }

    return result;
}

void parse_command_line(int argc, char* argv[], args_t& args) {

    try {
        for (int i = 1; i < argc; ++i) {
            if (std::string(argv[i]) == "--help") {
                std::cout << "Test script to compare serial and parallel merging. " << std::endl
                        << "You must add the path of the program 'benchmark_verify' " << std::endl
                        << "to your PATH environment variable." << std::endl << std::endl;
                std::cout << "Arguments:" << std::endl;
                std::cout << "\t"     << std::get<0>(args.nSamples)
                        << "\t\t"   << std::get<1>(args.nSamples) << std::endl;
                std::cout << "\t"     << std::get<0>(args.nRanks)
                        << "\t\t"   << std::get<1>(args.nRanks) << std::endl;
                std::cout << "\t"     << std::get<0>(args.nParcelPerCell)
                        << "\t\t\t" << std::get<1>(args.nParcelPerCell) << std::endl;
                std::cout << "\t"     << std::get<0>(args.minVratio)
                        << "\t\t"   << std::get<1>(args.minVratio) << std::endl;
                std::cout << "\t"     << std::get<0>(args.nx)
                        << "\t\t\t" << std::get<1>(args.nx) << std::endl;
                std::cout << "\t"     << std::get<0>(args.ny)
                        << "\t\t\t" << std::get<1>(args.ny) << std::endl;
                std::cout << "\t"     << std::get<0>(args.nz)
                        << "\t\t\t" << std::get<1>(args.ny) << std::endl;
                std::cout << "\t"     << std::get<0>(args.seed)
                        << "\t\t\t" << std::get<1>(args.seed) << std::endl;
                std::cout << "\t"     << std::get<0>(args.nTasksPerNode)
                        << "\t"     << std::get<1>(args.nTasksPerNode) << std::endl;
                std::cout << "\t"     << std::get<0>(args.commType)
                        << "\t\t"   << std::get<1>(args.commType) << std::endl;
                std::cout << "\t"     << std::get<0>(args.shuffle)
                        << "\t\t"   << std::get<1>(args.shuffle) << std::endl;
                std::cout << "\t"     << std::get<0>(args.subcomm)
                        << "\t\t"   << std::get<1>(args.subcomm) << std::endl;
                std::cout << "\t"     << std::get<0>(args.verbose)
                        << "\t\t"   << std::get<1>(args.verbose) << std::endl;
                std::cout << "\t"     << std::get<0>(args.cmd)
                        << "\t\t\t" << std::get<1>(args.cmd) << std::endl;
                return;
            } else if (std::string(argv[i]) == std::get<0>(args.nSamples)) {
                check_single_argument(argv, i);
                std::get<2>(args.nSamples) = std::atoi(argv[i+1]);
            } else if (std::string(argv[i]) == std::get<0>(args.nRanks)) {
                int nargs = check_multiple_argument(argv, i);
                std::list<int>& nRanks = std::get<2>(args.nRanks);
                for(int j = 1; j <= nargs; ++j) {
                    nRanks.push_back(std::atoi(argv[i+j]));
                }
            } else if (std::string(argv[i]) == std::get<0>(args.nParcelPerCell)) {
                check_single_argument(argv, i);
                std::get<2>(args.nParcelPerCell) = std::atoi(argv[i+1]);
            } else if (std::string(argv[i]) == std::get<0>(args.minVratio)) {
                check_single_argument(argv, i);
                std::get<2>(args.minVratio) = std::atof(argv[i+1]);
            } else if (std::string(argv[i]) == std::get<0>(args.nx)) {
                check_single_argument(argv, i);
                std::get<2>(args.nx) = std::atoi(argv[i+1]);
            } else if (std::string(argv[i]) == std::get<0>(args.ny)) {
                check_single_argument(argv, i);
                std::get<2>(args.ny) = std::atoi(argv[i+1]);
            } else if (std::string(argv[i]) == std::get<0>(args.nz)) {
                check_single_argument(argv, i);
                std::get<2>(args.nz) = std::atoi(argv[i+1]);
            } else if (std::string(argv[i]) == std::get<0>(args.seed)) {
                check_single_argument(argv, i);
                std::get<2>(args.seed) = std::atoi(argv[i+1]);
            } else if (std::string(argv[i]) == std::get<0>(args.nTasksPerNode)) {
                check_single_argument(argv, i);
                std::get<2>(args.nTasksPerNode) = std::atoi(argv[i+1]);
            } else if (std::string(argv[i]) == std::get<0>(args.commType)) {
                check_single_argument(argv, i);
                std::get<2>(args.commType) = std::string(argv[i+1]);
            } else if (std::string(argv[i]) == std::get<0>(args.shuffle)) {
                std::get<2>(args.shuffle) = true;
            } else if (std::string(argv[i]) == std::get<0>(args.verbose)) {
                std::get<2>(args.verbose) = true;
            } else if (std::string(argv[i]) == std::get<0>(args.cmd)) {
                check_single_argument(argv, i);
                std::get<2>(args.cmd) = std::string(argv[i+1]);
            }
        }
    } catch (std::exception& e) {
        std::cout << "Input error: " << e.what() << std::endl;
    }
}

int main(int argc, char* argv[]) {

    args_t args;

    parse_command_line(argc, argv, args);

    int nppc = std::get<2>(args.nParcelPerCell);
    int nx = std::get<2>(args.nx);
    int ny = std::get<2>(args.ny);
    int nz = std::get<2>(args.nz);
    int nParcels = nppc * nx * ny * nz;

    std::list<int> nRanks = std::get<2>(args.nRanks);

    double vcell = 1.0 / (nx * ny * nz);
    double minVratio = std::get<2>(args.minVratio);
    double vmin = vcell / minVratio;

    double tol = 1.0e-14;

    int nFails = 0;
    int nMerges = 0;
    int modulo = 100;

    int nTasksPerNode = std::get<2>(args.nTasksPerNode);
    int nSamples = std::get<2>(args.nSamples);
    int seed = std::get<2>(args.seed);
    bool verbose = std::get<2>(args.verbose);
    std::string commType = std::get<2>(args.commType);

    std::string exe;

    std::string flags = " --nx " + std::to_string(nx)
                      + " --ny " + std::to_string(ny)
                      + " --nz " + std::to_string(nz)
                      + " --n_per_cell " + std::to_string(nppc)
                      + " --min_vratio " + std::to_string(minVratio);

    if (std::get<2>(args.shuffle)) {
        flags = flags + " --shuffle";
    }

    std::string pflags = flags + " --comm-type " + commType;

    if (std::get<2>(args.subcomm)) {
        pflags = pflags + " --subcomm";
    }


    for (int n = 0; n < nSamples; ++n) {

        try {
            std::string cmd = "mpirun -np 1 ";
            if (std::get<2>(args.cmd) == "srun") {
                cmd = "srun --nodes=1 --ntasks=1 --ntasks-per-node=1 --exact ";
            }
            cmd = cmd + exe + " --seed " + std::to_string(seed) + " --setup-parcels";
            cmd = cmd + " > /dev/null 2>&1";
            std::system(cmd.c_str());
//             using namespace std::chrono_literals;
//             std::this_thread::sleep_for(120s);

        } catch(std::exception& e) {
            std::cout << "Error in running the serial version." << std::endl;
            return 0;
        }

        // Note: Because the 1 MPI run writes the initial parcel setup; the actual
        // solve has number 2 instead of 1.
        fs::path filename = "serial_final_0000000002_parcels.nc";
        if (fs::exists(filename)) {
            fs::rename(filename, "serial_final_0000000001_parcels.nc");

            std::cout << "Sample " << n << " generated." << std::endl;
        }

        for (int nRank : nRanks) {

            int nodes = int(nRank / nTasksPerNode + 0.5);

            // -------------------------------------------------------------
            // Run the serial and parallel versions of the nearest + merging algorithm:
            // We wait 2 minutes per process. If they exceed this time limit, we
            // assume the nearest algorithm is in an endless loop, i.e. deadlocked.
            bool failed = false;

            try {
                std::string cmd = "mpirun -np  " + std::to_string(nRank) + " ";
                if (std::get<2>(args.cmd) == "srun") {
                    cmd = "srun --nodes=" + std::to_string(nodes)
                        + " --ntasks=" + std::to_string(nRank);
                    cmd = cmd + " --cpus-per-task=1 --exact ";
                }
                cmd = cmd + exe + " " + pflags + " > /dev/null 2>&1";
                std::system(cmd.c_str());
//                 using namespace std::chrono_literals;
//                 std::this_thread::sleep_for(120s);

            } catch(std::exception& e) {
                std::cout << "Error in running the parallel version." << std::endl;
                failed = true;
            }

            failed = failed ||
                     !fs::exists("serial_final_0000000001_parcels.nc") ||
                     !fs::exists("parallel_final_0000000001_parcels.nc");

            if (!failed) {
                // ----------------------------------------
                // Compare the results:
                //FIXME

            }

            // -------------------------------------------------------------
            // Do clean up:
            if (failed) {
                failed = false;
                nFails = nFails + 1;
                std::string nstr = std::to_string(nFails);
                nstr.insert(0, 10-nstr.size(), '0');

                if (fs::exists("initial_0000000001_parcels.nc")) {
                    fs::rename("initial_0000000001_parcels.nc",
                               "initial_fail_" + nstr + "_parcels.nc");
                }

                if (fs::exists("serial_final_0000000001_parcels.nc")) {
                    fs::copy_file("serial_final_0000000001_parcels.nc",
                                  "serial_fail_" + nstr + "_parcels.nc");
                }

                if (fs::exists("parallel_final_0000000001_parcels.nc")) {
                    fs::rename("parallel_final_0000000001_parcels.nc",
                               "parallel_fail_" + nstr + "_parcels.nc");
                }
            } else {
                bool isRemoved = fs::remove("parallel_final_0000000001_parcels.nc");
                if (!isRemoved) {
                    std::cerr << "Unable to delete 'parallel_final_0000000001_parcels.nc'. "
                              << "File does not exist." << std::endl;
                }
            }
        }

        // -------------------------------------------------------------
        // Intermediate info:
        if (verbose && (n % modulo == 0)) {
            std::cout << "#samples, #fails, #merges: " << n << " " << nFails << " " << nMerges << std::endl;
        }

        bool isRemoved = fs::remove("initial_0000000001_parcels.nc");
        if (!isRemoved) {
            std::cerr << "Unable to delete 'initial_0000000001_parcels.nc'. "
                      << "File does not exist." << std::endl;
        }
        if (fs::exists("serial_final_0000000001_parcels.nc")) {
            isRemoved = fs::remove("serial_final_0000000001_parcels.nc");
            if (!isRemoved) {
                std::cerr << "Unable to delete 'serial_final_0000000001_parcels.nc'. "
                          << "File does not exist." << std::endl;
            }
        }

        seed++;
    }


    // -------------------------------------------------------------------------
    // Print summary:
    std::cout << "--------------------------------------------------------------------"
              << "Total number of samples:      " << nSamples << std::endl;
    std::cout << "MPI ranks:                    ";
    for (int nRank : nRanks) {
        std::cout << nRank << " ";
    }
    std::cout << std::endl
              << "Number of parcels per sample: " << nParcels << std::endl
              << "Number of parcels per cell:   " << nppc << std::endl
              << "Number of fails:              " << nFails << std::endl
              << "Number of merges:             " << nMerges << std::endl;

    return 0;
}
