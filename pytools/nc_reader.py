import netCDF4 as nc
from nc_base_reader import nc_base_reader
import re
import numpy as np
import os


class nc_reader(nc_base_reader):
    def __init__(self):
        # either a field or parcel file
        self._nctype = None

    def open(self, fname):
        super().open(fname)

        self._nctype = self.get_global_attribute("file_type")

        if not self._nctype == 'parcels':
            raise ValueError("No parcel file.")

        # if we read in a parcel file we pre-evaluate the number of
        # steps
        self._loaded_step = -1
        basename = os.path.basename(fname)
        # 14 Feb 2022
        # https://stackoverflow.com/questions/15340582/python-extract-pattern-matches
        p = re.compile(r"(.*)_(\d*)_parcels.nc")
        result = p.search(basename)
        self._basename = result.group(1)
        self._dirname = os.path.dirname(fname)
        if self._dirname == '':
            self._dirname = '.'
        self._n_parcel_files = 0
        for ff in os.listdir(self._dirname):
            if self._basename in ff and '_parcels.nc' in ff:
                self._n_parcel_files += 1

    def get_num_steps(self):
        return self._n_parcel_files

    def get_box_extent(self):
        return self.get_global_attribute("extent")

    def get_box_ncells(self):
        return self.get_global_attribute("ncells")

    def get_box_origin(self):
        return self.get_global_attribute("origin")

    def get_axis(self, name):
        if not name in ['x', 'y', 'z']:
            raise ValueError("No axis called '" + name + "'.")
        axis = self.get_all(name)
        if name == 'x' or name == 'y':
            axis = np.append(axis, abs(axis[0]))
        return axis

    def get_all(self, name):
        return super().get_all(name)

    def get_dataset(self, step, name, indices=None, copy_periodic=True):

        if name == 't':
            # parcel files store the time as a global attribute
            self._load_parcel_file(step)
            return np.array(self._ncfile.variables[name]).squeeze()

        self._load_parcel_file(step)
        if indices is not None:
            return np.array(self._ncfile.variables[name]).squeeze()[indices, ...]
        else:
            return np.array(self._ncfile.variables[name]).squeeze()

    def get_global_attribute(self, name):
        if not name in self._ncfile.ncattrs():
            raise IOError("Global attribute '" + name + "' unknown.")
        attr = self._ncfile.getncattr(name)
        if isinstance(attr, np.bytes_):
            attr = attr.decode()
        return attr

    def get_num_parcels(self, step):
        if not self.is_parcel_file:
            raise IOError("Not a parcel output file.")
        self._load_parcel_file(step)
        return self._ncfile.dimensions['n_parcels'].size

    def _get_step_string(self, step):
        return str(step).zfill(10)

    def _load_parcel_file(self, step):
        step = step + 1
        if self._loaded_step == step:
            return
        self._loaded_step = step
        self._ncfile.close()
        s = self._get_step_string(step)
        fname = os.path.join(self._dirname, self._basename + '_' + s + '_parcels.nc')
        self._ncfile = nc.Dataset(fname, "r", format="NETCDF4")
