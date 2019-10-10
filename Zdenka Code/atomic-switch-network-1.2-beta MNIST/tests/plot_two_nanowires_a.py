#!/usr/bin/python
"""
This demo plots the minimal unit of an atomic switch network: 
    + two wires and one junction.

This is a basic case in which the junction position coincides with the wires
centres.

.. moduleauthor:: Paula Sanz-Leon
August 2016 
"""
import numpy as np
from mayavi import mlab

mlab.figure(1, bgcolor=(0, 0, 0), size=(350, 350))
mlab.clf()

# Position of the extremes of the wires
extremes_x = np.array([1.0, 1.0, 0.0, 2.0]) * 40 / 5.5
extremes_y = np.array([0.0, 2.0, 1.0, 1.0]) * 40 / 5.5
extremes_z = np.array([0.0, 0.0, 1.0, 1.0]) * 40 / 5.5


# Position of the centres of the wires
centres_x = np.array([1.0, 1.0]) * 40 / 5.5
centres_y = np.array([1.0, 1.0]) * 40 / 5.5
centres_z = np.array([0.0, 1.0]) * 40 / 5.5

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

# The wires between two points in space -- this should be done with a dictionary
mlab.plot3d(extremes_x[0:2], extremes_y[0:2], extremes_z[0:2],
            tube_radius=0.5, color=(0.5, 0.5, 0.5))

mlab.plot3d(extremes_x[2:4], extremes_y[2:4], extremes_z[2:4], 
            tube_radius=0.5, color=(0.5, 0.5, 0.5))

# The junction
mlab.plot3d(centres_x, centres_y, centres_z, 
            tube_radius=0.4, color=(0.0, 0.0, 0.5), opacity=0.7)

mlab.show()