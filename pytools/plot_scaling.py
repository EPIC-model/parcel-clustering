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
    # Helper functions for parsing and plotting the data:
    #

    def get_configurations(path, test_case):
        tc = '-' + test_case + '-'
        pattern = re.compile(r"(\w*)-(\w*)" + tc + r"(nx-\d*-ny-\d*)-nodes-(\d*).nc")

        compilers = set()
        grids = set()
        graphs = set()

        groups = set()
        configs = {}
        for fname in os.listdir(path=path):
            m = re.match(pattern, fname)
            if not m is None:
                group = m.group(1) + '-' + m.group(2) + tc + m.group(3)
                groups.add(group)
                compilers.add(m.group(1))
                graphs.add(m.group(2))
                grids.add(m.group(3))
                if not group in configs.keys():
                    configs[group] = {
                        'basename': group + '-nodes-',
                        'compiler': m.group(1),
                        'graph':    m.group(2),
                        'grid':     m.group(3),
                        'nodes':    []
                    }
                configs[group]['nodes'].append(int(m.group(4)))

        print("Found", len(configs.keys()), "different configurations. There are")
        print("\t", len(compilers), "compiler(s):", compilers)
        print("\t", len(graphs), "graph method(s):", graphs)
        print("\t", len(grids), "grid configuration(s):", grids)
        print()

        # sort nodes
        for group in groups:
            if not group in configs.keys():
                raise KeyError("Group '" + group + "' not found.")
            config = configs[group]
            config['nodes'].sort()

        return configs

    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    def add_to_plot(ax, config, timings, cmap, marker, add_label=False):
        nodes = config['nodes']

        # switch case:
        # (see https://docs.python.org/3.10/whatsnew/3.10.html#pep-634-structural-pattern-matching, 28 Jan 2025)
        method = ''
        match config['graph']:
            case 'caf':
                method = 'CAF'
            case 'p2p':
                method = 'MPI-3 P2P'
            case 'rma':
                method = 'MPI-3 RMA'
            case 'shmem':
                method = 'SHMEM'
            case _:
                raise RuntimeError("No method called '" + config['graph'] + "'.")
        # done

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

    def generate_plot(configs, comm, args):

        tag = args.compiler_suite + '-' + comm + '-' + args.test_case

        titles = {
            'p2p':   'MPI point-to-point communication',
            'rma':   'MPI RMA communication',
            'shmem': 'SHMEM communication',
            'caf':   'Coarray Fortran'
        }

        title = titles[comm]

        cmap = plt.get_cmap(args.colour_map)

        markers = args.markers

        groups = list(configs.keys())

        n_conf = sum(tag in group for group in groups)

        if n_conf > len(markers):
            raise RuntimeError('Not enough markers. ' + \
                'Please add more to the command line with --markers')

        plt.figure(figsize=(8, 7), dpi=200)

        ax = plt.gca()


        ax.set_title(title)
        ax.grid(which='both', linestyle='dashed', linewidth=0.25)


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

        # -----------------------------------------------------------
        # Create legend where markers share a single legend entry:
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
                  handler_map={list: HandlerTuple(ndivide=None)})

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
        default='all',
        choices=['all', 'p2p', 'rma', 'shmem', 'caf'],
        help="Communication method.",
    )

    parser.add_argument(
        "--timings",
        type=str,
        nargs='+',
        default=['parcel merge (total)', 'find nearest', 'resolve graphs'],
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

    args = parser.parse_args()

    configs = get_configurations(args.path, args.test_case)

    if args.comm == 'all':
        for comm in ['p2p', 'rma', 'shmem', 'caf']:
            generate_plot(configs, comm=comm, args=args)
    else:
        generate_plot(configs, comm=args.comm, args=args)

except Exception as ex:
    print(ex, flush=True)

