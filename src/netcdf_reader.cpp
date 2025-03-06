#include "netcdf_reader.h"
#include <exception>
#include <iostream>

NetcdfReader::NetcdfReader() : ncid_m(-1)
{ }

bool NetcdfReader::open(std::filesystem::path filename) {
    int status = nc_open(filename.string().c_str(), NC_NOWRITE, &this->ncid_m);
    return (status == NC_NOERR);
}

bool NetcdfReader::close() {
    int status = nc_close(this->ncid_m);
    return (status == NC_NOERR);
}

size_t NetcdfReader::get_dimension(std::string var) const {
    int varid = -1;
    int status = nc_inq_dimid(this->ncid_m, var.c_str(), &varid);

    if (status != NC_NOERR) {
        std::runtime_error("NetCDF: Unable to find '" + var + "'.");
    }

    size_t varlength = 0;
    status = nc_inq_dimlen(this->ncid_m, varid, &varlength);

    if (status != NC_NOERR) {
        std::runtime_error("NetCDF: Unable to retrieve dimension of '" + var + "'.");
    }

    return varlength;
}

std::vector<double> NetcdfReader::get_dataset(std::string var) const {

    int varid = -1;
    int status = nc_inq_varid(this->ncid_m, var.c_str(), &varid);

    if (status != NC_NOERR) {
        std::runtime_error("NetCDF: Unable to find dataset '" + var + "'.");
    }

    size_t count = this->get_dimension("n_parcels");

    std::vector<double> values(count);

    status = nc_get_var_double(this->ncid_m, varid, &values[0]);

    if (status != NC_NOERR) {
        std::runtime_error("NetCDF: Unable to read data of '" + var + "'.");
    }

    return values;
}
