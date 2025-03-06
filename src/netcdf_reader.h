#ifndef NETCDF_READER
#define NETCDF_READER

#include <filesystem>
#include <netcdf.h>
#include <vector>

class NetcdfReader {

public:

    NetcdfReader();

    bool open(std::filesystem::path filename);

    bool close();


    size_t get_dimension(std::string var) const;

    std::vector<double> get_dataset(std::string var) const;

private:
    int ncid_m;

};

#endif
