// g++ -std=c++23 -I $NETCDF_C_INCLUDE_DIR netcdf_reader.h netcdf_reader.cpp verify_cluster_algorithm.cpp -L $NETCDF_C_LIBRARY_DIR -lnetcdf
#include <iostream>
#include <string>
#include <tuple>
#include <list>
#include <cstdlib>
#include <filesystem>
#include <exception>
#include <numeric>
#include <algorithm>
#include <vector>
#include <cmath>

#include "netcdf_reader.h"

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
    if (argv[i+1] == nullptr || std::string(argv[i+1]).find("--") != std::string::npos) {
        throw std::invalid_argument("No value provided for " + std::string(argv[i]));
    }
}

int check_multiple_argument(char* argv[], int i) {
    int result = 0;

    int j = i+1;
    while (argv[j] != nullptr && !std::string(argv[j]).find("--") != std::string::npos) {
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

void sort(const NetcdfReader& ncr, std::vector<int>& indices) {

    std::vector<double> x = ncr.get_dataset("x_position");
    std::vector<double> y = ncr.get_dataset("y_position");
    std::vector<double> z = ncr.get_dataset("z_position");

    indices.clear();

    indices.resize(x.size());
    std::iota(indices.begin(), indices.end(), 0);

    std::stable_sort(indices.begin(), indices.end(),
                     [&x, &y, &z](int i, int j) {
                         return (x[i] < x[j]) && (y[i] < y[j]) && (z[i] < z[j]);
                    });
}


bool compare_results(int nRank) {
    double tol = 1.0e-14;

    NetcdfReader ncrs, ncrp;

    if (!ncrs.open("serial_final_0000000001_parcels.nc")) {
        throw std::runtime_error("Unable to open 'serial_final_0000000001_parcels.nc'.");
    }

    if (!ncrp.open("parallel_final_0000000001_parcels.nc")) {
        throw std::runtime_error("Unable to open 'parallel_final_0000000001_parcels.nc'.");
    }


    size_t world_size = ncrp.get_dimension("world%size");

    if (nRank != world_size) {
        std::cout << "Warning: Failed to run with requested number of MPI cores"
                  << std::endl << std::flush;
    }

    size_t count = ncrs.get_dimension("n_parcels");

    if (count != ncrp.get_dimension("n_parcels")) {
        return true;
    }

    std::vector<int> ind1, ind2;
    sort(ncrs, ind1);
    sort(ncrp, ind2);

    std::list<std::string> attributes =
    {
        "x_position", "y_position", "z_position",
        "x_vorticity", "y_vorticity", "z_vorticity",
        "buoyancy", "volume", "B11", "B12", "B13", "B22", "B23"
    };

    for (const std::string& attr : attributes) {

        std::vector<double> ds1 = ncrs.get_dataset(attr);
        std::vector<double> ds2 = ncrp.get_dataset(attr);

        double error = 0.0;

        for(size_t j = 0; j < count; ++j) {
            int k = ind1[j];
            int l = ind2[j];
            error = std::max(error, std::abs(ds1[k] - ds2[l]));
        }

        if (error > tol) {
            return true;
        }
    }


    if (!ncrs.close()) {
        throw std::runtime_error("Unable to close 'serial_final_0000000001_parcels.nc'.");
    }

    if (!ncrp.close()) {
        throw std::runtime_error("Unable to close 'serial_final_0000000001_parcels.nc'.");
    }

    return false;
}

int main(int argc, char* argv[]) {

    try {
        args_t args;

        parse_command_line(argc, argv, args);

        const int nppc = std::get<2>(args.nParcelPerCell);
        const int nx = std::get<2>(args.nx);
        const int ny = std::get<2>(args.ny);
        const int nz = std::get<2>(args.nz);
        const int nParcels = nppc * nx * ny * nz;

        const std::list<int> nRanks = std::get<2>(args.nRanks);

        const double vcell = 1.0 / (nx * ny * nz);
        const double minVratio = std::get<2>(args.minVratio);
        const double vmin = vcell / minVratio;

        int nFails = 0;
        int nMerges = 0;
        int modulo = 100;

        const int nTasksPerNode = std::get<2>(args.nTasksPerNode);
        const int nSamples = std::get<2>(args.nSamples);
        int seed = std::get<2>(args.seed);
        const bool verbose = std::get<2>(args.verbose);
        const std::string commType = std::get<2>(args.commType);

        const std::string exe = "benchmark_verify";

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

            } catch(const std::exception& e) {
                throw std::runtime_error("Error in running the serial version.");
            }

            if (!fs:exists("initial_0000000001_parcels.nc")) {
                throw std::runtime_error("Input netCDF file not generated.");
            }

            // Note: Because the 1 MPI run writes the initial parcel setup; the actual
            // solve has number 2 instead of 1.
            fs::path filename = "serial_final_0000000002_parcels.nc";
            if (!fs::exists(filename)) {
                throw std::runtime_error("Serial output netCDF file not generated.");
            }

            fs::rename(filename, "serial_final_0000000001_parcels.nc");

            std::cout << "Sample " << n << " generated." << std::endl << std::flush;

            for (int nRank : nRanks) {

                const int nodes = int(std::ceil(double(nRank) / double(nTasksPerNode)));

                try {
                    std::string cmd = "mpirun -np  " + std::to_string(nRank) + " ";
                    if (std::get<2>(args.cmd) == "srun") {
                        cmd = "srun --nodes=" + std::to_string(nodes)
                            + " --ntasks=" + std::to_string(nRank);
                        cmd = cmd + " --cpus-per-task=1 --exact ";
                    }
                    cmd = cmd + exe + " " + pflags + " > /dev/null 2>&1";
                    std::system(cmd.c_str());

                } catch(const std::exception& e) {
                    throw std::runtime_error("Error in running the parallel version.");
                }

                if (!fs::exists("parallel_final_0000000001_parcels.nc")) {
                    throw std::runtime_error("Parallel output netCDF file not generated.");
                }

                bool failed = compare_results(nRank);

                // -------------------------------------------------------------
                // Do clean up:
                if (failed) {
                    nFails = nFails + 1;
                    std::string nstr = std::to_string(nFails);
                    nstr.insert(0, 10-nstr.size(), '0');

                    fs::copy_file("initial_0000000001_parcels.nc",
                                  "initial_fail_" + nstr + "_parcels.nc");

                    fs::copy_file("serial_final_0000000001_parcels.nc",
                                  "serial_fail_" + nstr + "_parcels.nc");

                    fs::rename("parallel_final_0000000001_parcels.nc",
                               "parallel_" + std::to_string(nRank) + "_fail_" + nstr + "_parcels.nc");

                } else {
                    bool isRemoved = fs::remove("parallel_final_0000000001_parcels.nc");
                    if (!isRemoved) {
                        throw std::runtime_error(
                            "Unable to delete 'parallel_final_0000000001_parcels.nc'.");
                    }
                }
            }

            // -------------------------------------------------------------
            // Intermediate info:
            if (verbose && (n % modulo == 0)) {
                std::cout << "#samples, #fails, #merges: "
                          << n << " " << nFails << " " << nMerges
                          << std::endl << std::flush;
            }

            bool isRemoved = fs::remove("initial_0000000001_parcels.nc");
            if (!isRemoved) {
                throw std::runtime_error(
                    "Unable to delete 'initial_0000000001_parcels.nc'.");
            }

            isRemoved = fs::remove("serial_final_0000000001_parcels.nc");
            if (!isRemoved) {
                throw std::runtime_error(
                    "Unable to delete 'serial_final_0000000001_parcels.nc'.");
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

    } catch(const std::exception& e) {
        std::cout << e.what() << std::endl;
    }

    return 0;
}
