#!/usr/bin/python
"""
This demo plots the minimal unit of an atomic switch network: 
    + two wires and one junction.

The junction position does not match the wires
centres.

.. moduleauthor:: Paula Sanz-Leon
August 2016 
"""
import sys
import numpy as np
from mayavi import mlab

# Append path
sys.path.append('/home/paula/Work/Code/AtomicSwitchNetworks/atomic-switch-network')
# Load functions to find intersection between two lines
from connectivity import find_intersection as fint


mlab.figure(1, bgcolor=(0, 0, 0), size=(350, 350))
mlab.clf()

# Position of the extremes of the wires
extremes_x = np.array([1.0, 1.0, 0.5, 1.5]) * 40 / 5.5
extremes_y = np.array([0.0, 3.0, 1.5, 2.5]) * 40 / 5.5
extremes_z = np.array([0.0, 0.0, 1.0, 1.0]) * 40 / 5.5

# Position of the centres of the wires
centres_x = np.array([1.0, 1.0]) * 40 / 5.5
centres_y = np.array([1.5, 2.0]) * 40 / 5.5
centres_z = np.array([0.0, 1.0]) * 40 / 5.5

# Find intersection between wires
# W = line([xa,ya], [xb, yb])
W1 = fint.line([extremes_x[0],extremes_y[0]], [extremes_x[1],extremes_y[1]])
W2 = fint.line([extremes_x[2],extremes_y[2]], [extremes_x[3],extremes_y[3]])

# Get xy coordinates of the junction
J_xy  = np.array(fint.intersection(W1, W2))

# Repeat junction xy coordinates 
junction_x = np.repeat(np.atleast_2d(J_xy), 2, axis=0)[:, 0]
junction_y = np.repeat(np.atleast_2d(J_xy), 2, axis=0)[:, 1]
junction_z = np.array([0.0, 1.0]) * 40 / 5.5

# Plot balls at the end points
AB = mlab.points3d(extremes_x, extremes_y, extremes_z,
                  scale_factor=1,
                  resolution=20,
                  color=(0.5, 0.5, 0.5),
                  scale_mode='none')

# Plot balls to mark the centres
C = mlab.points3d(centres_x, centres_y, centres_z,
                  scale_factor=1,
                  resolution=20,
                  color=(1, 0, 0),
                  opacity=0.5,
                  scale_mode='none')

# Plot balls to mark the junction
J = mlab.points3d(junction_x, junction_y, junction_z,
                  scale_factor=1,
                  resolution=20,
                  color=(0, 0, 1),
                  opacity=0.7,
                  scale_mode='none')

# The wires between two points in space -- this should be done with a dictionary
mlab.plot3d(extremes_x[0:2], extremes_y[0:2], extremes_z[0:2],
            tube_radius=0.5, color=(0.5, 0.5, 0.5))

mlab.plot3d(extremes_x[2:4], extremes_y[2:4], extremes_z[2:4], 
            tube_radius=0.5, color=(0.5, 0.5, 0.5))

# The junction
mlab.plot3d(junction_x, junction_y, junction_z, 
            tube_radius=0.4, color=(0.0, 0.0, 0.5), opacity=0.7)

mlab.show()