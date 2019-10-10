#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This demo shows how to produce a network of nanowires using 
the functions in the module `wires`.

This script should be used to produce networks with different
properties.

Assumes you have the module wires.py in your PYTHONPATH

For help type:
    python generate_nanowires_network.py -h

OPTIONS:
    python generate_nanowires_network.py 

or,
    python  generate_nanowires_network.py --nwires 100 --mean_length 100.0 --std_length 10.0 --seed 77

.. moduleauthor:: Paula Sanz-Leon <paula.sanz-leon@sydney.edu.au>
"""    

import wires

nwires = 995
centroid_dispersion = 700
mean_length = 100.0
std_length = 10.0
seed = 42
Lx   = 3000
Ly   = 3000

plot_network = 0

import logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s:%(levelname)s:%(message)s')


# Generate the network
wires_dict = wires.generate_wires_distribution(number_of_wires = nwires,
                                         wire_av_length = mean_length,
                                         wire_dispersion = std_length,
                                         gennorm_shape = 3,
                                         centroid_dispersion=centroid_dispersion,
                                         this_seed = seed,
                                         Lx = Lx,
                                         Ly = Ly)

# Get junctions list and their positions
wires.detect_junctions(wires_dict)

# Genreate graph object and adjacency matrix
wires.generate_graph(wires_dict)

if not wires.check_connectedness(wires_dict):
    logging.warning("Will select the largest connected component.")
    wires_dict = wires.select_largest_component(wires_dict)

    logging.info("The graph is connected. Will save it to mat file.")
    wires.export_to_matlab(wires_dict)

if plot_network:
 
    # Plotting tools
    from matplotlib.lines import Line2D
    from matplotlib.patches import Rectangle
    import matplotlib.pyplot as plt

    # Plot pretty pictures of what we just did
    fig, ax = plt.subplots()
    fig.set_size_inches(5,5)

    Lx = wires_dict['length_x']
    Ly = wires_dict['length_y']

    ax.add_patch(Rectangle((0,0), Lx, Ly, color=(1.0, 0.918, 0.0), alpha=0.77))     
    ax = draw_wires(ax, wires_dict)
    ax = draw_junctions(ax, wires_dict)
    ax.set_aspect(1) # set aspect ratio to 1
    ax.set_xlabel(r'x [$\mu$m]')
    ax.set_ylabel(r'y [$\mu$m]')
    ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    ax.axis([-.1*Lx,1.1*Lx,-.1*Lx,1.1*Lx]) # add some space around the unit square
    ax.set_title('Nanowires distribution')
    ax.grid()
    plt.show()

