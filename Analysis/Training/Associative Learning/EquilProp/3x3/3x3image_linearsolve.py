
# -*- coding: utf-8 -*-
"""
Created on Thu May 27 00:36:08 2021

@author: gnate
"""

from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from numpy import linalg
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.lines as mlines
import pandas as pd
import skspatial
from skspatial.objects import Plane, Points
from skspatial.plotting import plot_3d

# Define function to generate one-hot encoding
def create_target(t):
    target_vector=np.zeros([t.shape[0],3])
    for row in range(t.shape[0]):
        target_vector[row,:] = np.zeros(3)
        for i in range(3):
            if i == t[row]:
                   target_vector[row,i] = 1    
    return target_vector

def create_single_target(t):
        target_vector = np.zeros(3)
        for i in range(3):
            if i == t:
                target_vector[i] = 1
        return target_vector
    
# Define angle function for rotating plot gif    
# def rotate(angle):
#        ax.view_init(azim=angle)


# Load 3x3 dataset
letters = pd.read_csv("./3x3image_gen.csv", sep=',', header = None)
data = letters.values[:, 0:9]

sample_num = 30      # Number of data samples

# View a random data sample as 3x3 image
# rnd = np.random.randint(0, 30)
# image_example = data[rnd].reshape(3,3).T

# Target labels: 0 - 'z', 1 - 'v', 2 - 'n'
targets = letters.values[:, 9]
inputs = np.zeros((sample_num,10))
inputs[:,0:9] = data
inputs[:,9] = np.ones(sample_num)   #Bias input

# Standardize data - Uncomment to standardize data, not necessary
# inputs = data - np.mean(data)
# inputs = inputs/(np.std(data))

# One-hot encoding of outputs
onehot_outputs = create_target(targets)

# Generate inverse of input matrix
inv_inputs = np.linalg.pinv(inputs)

# Linear solve for inputs * W = onehot_outputs
W = np.zeros((10,3))
W = np.matmul(inv_inputs, onehot_outputs)
print("Weight matrix from linear solve - \n", W[:9,:])
print("Biases of output nodes - \n", W[9,:])

# Generate scatter plot of predicted outputs from linear solve 
predicted_outputs = np.matmul(inputs, W)
#print("\n Output predictions using linear solve - \n", predicted_outputs)

# Calculate prediction accuracy
accuracy = 0
for i in range(0,sample_num):
    if np.argmax(predicted_outputs[i]) == targets[i]:
        accuracy = accuracy + 1
print('Prediction accuracy -', accuracy/sample_num)

# Plot figures 
fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111, projection='3d')

for c, m, low, high in [('r', 'o', 0, 10), ('b', '^', 10, 20), ('g', 'v', 20, sample_num)]:
    xdata = predicted_outputs[low:high,0]
    ydata = predicted_outputs[low:high,1]
    zdata = predicted_outputs[low:high,2]
    ax.scatter(xdata, ydata, zdata, c = c, marker = m);
    
# Make simple, bare axis lines through space:
xAxisLine = ((min(predicted_outputs[:,0]), max(predicted_outputs[:,0])), (0, 0), (0,0))
ax.plot(xAxisLine[0], xAxisLine[1], xAxisLine[2], 'r')
yAxisLine = ((0, 0), (min(predicted_outputs[:,1]), max(predicted_outputs[:,1])), (0,0))
ax.plot(yAxisLine[0], yAxisLine[1], yAxisLine[2], 'r')
zAxisLine = ((0, 0), (0,0), (min(predicted_outputs[:,2]), max(predicted_outputs[:,2])))
ax.plot(zAxisLine[0], zAxisLine[1], zAxisLine[2], 'r')

# Set axes labels and figure legend
ax.set_xlabel('X-axis')
ax.set_ylabel('Y-axis')
ax.set_zlabel('Z-axis')

red_point =  mlines.Line2D([], [], color='red', marker='o', 
                          markersize=10, label='Z')
blue_point = mlines.Line2D([], [], color='blue', marker='^', 
                           markersize=10, label='V')
green_point = mlines.Line2D([], [], color='green', marker='v', 
                           markersize=10, label='N')
plt.legend(handles=[red_point, blue_point, green_point])

plt.show()

# Generate rotating gif - Uncomment to generate the gif in current dir
# angle = 3
# ani = animation.FuncAnimation(fig, rotate, frames=np.arange(0, 360, angle), interval=20)
# ani.save('predicted_outputs.gif', writer=animation.PillowWriter(fps=10))

# Test with some input
#rnd = np.random.randint(0,29)
#print("Predicted output -", test_output)
#actual_output = create_single_target(targets[rnd])
#print("Actual output -", actual_output)



