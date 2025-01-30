import netCDF4 as nc
import numpy as np
import os
import matplotlib.pyplot as plt
import re
import argparse
from matplotlib.legend_handler import HandlerTuple

plt.rcParams['font.family'] = 'sans'
plt.rcParams['font.size'] = 12
plt.rcParams['text.usetex'] = True

try:

    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # Helper functions and classes for parsing and plotting the data:
    #
    class DataSet:

        def __init__(self, path, test_case):
            self.compilers = set()
            self.grids = set()
            self.comms = set()
            self.groups = set()

            self.titles = {
                'p2p':   r'MPI-3 P2P', # communication',
                'rma':   r'MPI-3 RMA', # communication',
                'shmem': r'SHMEM', # communication',
                'caf':   r'Coarray Fortran'
            }

            tc = '-' + test_case + '-'
            pattern = re.compile(r"(\w*)-(\w*)" + tc + r"(nx-\d*-ny-\d*)-nodes-(\d*).nc")

            self.configs = {}
            for fname in os.listdir(path=path):
                m = re.match(pattern, fname)
                if not m is None:
                    group = m.group(1) + '-' + m.group(2) + tc + m.group(3)
                    self.groups.add(group)
                    self.compilers.add(m.group(1))
                    self.comms.add(m.group(2))
                    self.grids.add(m.group(3))
                    if not group in self.configs.keys():
                        self.configs[group] = {
                            'basename': group + '-nodes-',
                            'compiler': m.group(1),
                            'comm':     m.group(2),
                            'grid':     m.group(3),
                            'nodes':    []
                        }
                    self.configs[group]['nodes'].append(int(m.group(4)))

            print("Found", len(self.configs.keys()), "different configurations. There are")
            print("\t", len(self.compilers), "compiler(s):", self.compilers)
            print("\t", len(self.comms), "comm method(s):", self.comms)
            print("\t", len(self.grids), "grid configuration(s):", self.grids)
            print()

            # sort nodes
            for group in self.groups:
                if not group in self.configs.keys():
                    raise KeyError("Group '" + group + "' not found.")
                config = self.configs[group]
                config['nodes'].sort()


        def get_data(self, config, nodes, timings):
            avg_data = {}
            std_data = {}
            for long_name in timings:
                avg_data[long_name] = np.zeros(len(nodes))
                std_data[long_name] = np.zeros(len(nodes))

            for i, node in enumerate(nodes):
                fname = config['basename'] + str(node) + '.nc'

                nc_file = nc.Dataset(fname, "r", format="NETCDF4")

                lvars = list(nc_file.variables.keys())

                for var in lvars:
                    name = nc_file.variables[var].name
                    if 'wtime' in name:
                        long_name = nc_file.variables[var].long_name

                        if long_name in timings:
                            data = np.array(nc_file[name])
                            avg_data[long_name][i] = data.mean()
                            std_data[long_name][i] = data.std()

                nc_file.close()
            return avg_data, std_data

    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    def add_to_plot(ax, config, timings, cmap, marker, add_label=False):
        nodes = np.asarray(config['nodes'])

        # switch case:
        # (see https://docs.python.org/3.10/whatsnew/3.10.html#pep-634-structural-pattern-matching, 28 Jan 2025)
        method = ''
        match config['comm']:
            case 'caf':
                method = 'CAF'
            case 'p2p':
                method = 'MPI-3 P2P'
            case 'rma':
                method = 'MPI-3 RMA'
            case 'shmem':
                method = 'SHMEM'
            case _:
                raise RuntimeError("No method called '" + config['comm'] + "'.")
        # done

        avg_data, std_data = dset.get_data(config, nodes, timings)

        label = None
        if add_label:
            label = 'ideal scaling'

        ax.plot(nodes,
                avg_data[timings[0]][0] / nodes * nodes[0],
                color='black',
                linestyle='dashed',
                linewidth=1,
                label=label)
        for i, long_name in enumerate(timings):
            label = None
            if add_label:
                label = long_name
            if "resolve graphs" in label:
                label = label + ' (' + method + ')'
            ax.errorbar(x=nodes,
                        y=avg_data[long_name],
                        yerr=std_data[long_name],
                        #yerr=abs(avg_data[long_name]-std_data[long_name]),
                        label=label,
                        color=cmap(i),
                        linewidth=1,
                        marker=marker,
                        markersize=5,
                        capsize=3)

    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    # Create legend where markers share a single legend entry:
    def add_legend(ax, **kwargs):
        handles, labels = ax.get_legend_handles_labels()

        arg = {}
        arg['ideal scaling'] = None
        for t in labels:
            arg[t] = []

        for i in range(len(labels)):
            if labels[i] == 'ideal scaling':
                arg['ideal scaling'] = handles[i]
            else:
                arg[labels[i]].append(handles[i])

        h_ = []
        for l in arg.keys():
            h_.append(arg[l])

        # 23 Jan 2025
        # https://matplotlib.org/stable/gallery/text_labels_and_annotations/legend_demo.html
        ax.legend(loc='lower left', handles=h_, labels=arg.keys(),
                  handler_map={list: HandlerTuple(ndivide=None)},
                  **kwargs)

    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    def generate_scaling_plot(dset, args):

        if 'all' in args.comm:
            n = len(dset.comms)
            nrows = int(np.sqrt(n))
            ncols = int(n / nrows + 0.5)

            sharex = (nrows > 1)
            fig, axs = plt.subplots(nrows=nrows,
                                    ncols=ncols,
                                    sharey=True,
                                    sharex=sharex,
                                    figsize=(5*ncols, 5*nrows),
                                    dpi=400)

            axs = axs.flatten()
            comms = sorted(dset.comms)
            for i, comm in enumerate(comms):
                axs[i].grid(which='both', linestyle='dashed', linewidth=0.25)

                # -----------------------------------------------------------
                # Add individual scaling:
                tag = args.compiler_suite + '-' + comm + '-' + args.test_case
                add_line(axs[i], dset, tag, comm, args)

                # -----------------------------------------------------------
                # Create legend where markers share a single legend entry:
                add_legend(axs[i],
                           title=r'\bfseries{' + dset.titles[comm] + r'}',
                           alignment='left')

                axs[i].set_yscale('log', base=10)
                axs[i].set_xscale('log', base=2)

                if i >= (nrows - 1) * ncols:
                    axs[i].set_xlabel('number of nodes (1 node = 128 cores)')

            axs[0].set_ylabel('run time (s)')

            # -----------------------------------------------------------
            # Save figure:
            plt.tight_layout()
            plt.savefig(args.compiler_suite + '-scaling.pdf', bbox_inches='tight')
            plt.close()

        else:
            for comm in args.comm:
                # -----------------------------------------------------------
                # Create figure:
                plt.figure(figsize=(8, 7), dpi=200)
                ax = plt.gca()
                title = dset.titles[comm]
                ax.set_title(title)
                ax.grid(which='both', linestyle='dashed', linewidth=0.25)

                # -----------------------------------------------------------
                # Add individual scaling:
                tag = args.compiler_suite + '-' + comm + '-' + args.test_case

                add_line(ax, dset, tag, comm, args)

                # -----------------------------------------------------------
                # Create legend where markers share a single legend entry:
                add_legend(ax)

                ax.set_yscale('log', base=10)
                ax.set_xscale('log', base=2)
                ax.set_xlabel('number of nodes (1 node = 128 cores)')
                ax.set_ylabel('run time (s)')

                # -----------------------------------------------------------
                # Save figure:
                plt.tight_layout()
                plt.savefig(tag + '.pdf', bbox_inches='tight')
                plt.close()

    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    def generate_efficiency_plot(dset, args):

        if 'all' in args.comm:
            raise RuntimeError("Not yet implemented.")
        else:
            for grid in sorted(dset.grids):
                # -----------------------------------------------------------
                # Create figure:
                plt.figure(figsize=(8, 7), dpi=200)
                ax = plt.gca()
                ax.grid(which='both', linestyle='dashed', linewidth=0.25, axis='y')

                # -----------------------------------------------------------
                # Add individual scaling:
                comms = sorted(args.comm)
                n_comms = len(comms)
                width= 0.4 / n_comms
                for i, comm in enumerate(comms):
                    tag = args.compiler_suite + '-' + comm + '-' + args.test_case + '-' + grid

                    offset = width * (i - 0.5*n_comms)
                    add_bar(ax, dset, tag, comm, args, offset=offset, width=width)


                ax.legend(loc='upper left')

                ax.axhline(y=1, linestyle='dashed', color='black')

                ax.set_xlabel('number of nodes (1 node = 128 cores)')
                ax.set_ylabel('strong efficiency')

                # -----------------------------------------------------------
                # Save figure:
                plt.tight_layout()
                plt.savefig(args.compiler_suite + '-' + grid + '-strong-effiency.pdf',
                            bbox_inches='tight')
                plt.close()

    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    def add_line(ax, dset, tag, comm, args):

        configs = dset.configs

        cmap = plt.get_cmap(args.colour_map)

        markers = args.markers

        groups = list(configs.keys())

        n_conf = sum(tag in group for group in groups)

        if n_conf > len(markers):
            raise RuntimeError('Not enough markers. ' + \
                'Please add more to the command line with --markers')

        found = True
        i = 0
        for group in groups:

            if not tag in group:
                continue

            config = configs[group]

            add_to_plot(ax,
                        config,
                        timings=args.timings,
                        cmap=cmap,
                        marker=markers[i],
                        add_label=True)

            i = i + 1
            found = False

    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    def add_bar(ax, dset, tag, comm, args, offset, width):

        configs = dset.configs

        cmap = plt.get_cmap(args.colour_map)

        groups = list(configs.keys())

        label = dset.titles[comm]

        for group in groups:

            if not tag in group:
                continue

            config = configs[group]

            nodes = np.asarray(config['nodes'])

            avg_data, std_data = dset.get_data(config, nodes, args.timings)

            speedup = avg_data['parcel merge'][0] / avg_data['parcel merge']
            p = nodes / nodes[0]
            eff = speedup / p
            n = len(p)
            x = np.arange(n)
            ax.bar(x=x+offset,
                   height=eff,
                   width=width,
                   align='edge',
                   label=label)

            ax.set_xticks(x, nodes)

    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # Actual 'main':
    #

    parser = argparse.ArgumentParser(
            description="Generate strong and weak scaling plot."
    )

    parser.add_argument(
        "--compiler-suite",
        type=str,
        default="cray",
        choices=['cray', 'gnu'],
        help="Compiler environment",
    )

    parser.add_argument(
        "--test-case",
        type=str,
        default="random",
        choices=['random', 'read'],
        help="Test case to analyse.",
    )

    parser.add_argument(
        "--comm",
        type=str,
        nargs='+',
        default='all',
        choices=['all', 'p2p', 'rma', 'shmem', 'caf'],
        help="Communication method.",
    )

    parser.add_argument(
        "--timings",
        type=str,
        nargs='+',
        default=['parcel merge', 'merge nearest', 'graph resolve'],
        #default=['parcel merge (total)', 'find nearest', 'build graphs', 'resolve graphs'],
        help="Timer data to visualise.",
    )

    parser.add_argument(
        "--path",
        type=str,
        default='.',
        help="Data directory.",
    )

    parser.add_argument(
        "--colour-map",
        type=str,
        default='tab10',
        help="Colour map for plotting."
    )

    parser.add_argument(
        "--markers",
        type=str,
        nargs='+',
        default=['o', 's', 'D'],
        help="Markers for line plot."
    )

    parser.add_argument(
        "--plot",
        type=str,
        default='scaling',
        choices=['scaling', 'weak-efficiency', 'strong-efficiency'],
        help="Plot scaling or efficiency figures.")

    args = parser.parse_args()

    dset = DataSet(args.path, args.test_case)

    if args.plot == 'scaling':
        generate_scaling_plot(dset, args=args)
    else:
        generate_efficiency_plot(dset, args=args)

except Exception as ex:
    print(ex, flush=True)

