import analyse_network
import scipy.io
import numpy as np
import os

os.chdir('connectivity_data\\network_data-23-2-2017') #will only work on windows, not sure

params_dict = {} #empty dict where we will store network params
print 'there are ' + str(len(os.listdir('.'))) + ' networks to analyse'

#loop through each file name
for i in os.listdir('.'):
    #load file and tie each set of parameters to their seed
    wires_dict = scipy.io.loadmat(i)
    seed = wires_dict['this_seed']
    params_dict[int(seed)] = analyse_network.caculate_graph_params(i)

#initialise dictionaries for storing the average and standard devaition of network params
avgs_dict = {}
std_dict = {}

# assumes that all seeds calculated same network params, takes keys from first seed dict and itterates over every key in that seed

for i in params_dict[params_dict.keys()[0]]: 
    avgs_dict[i] = 0.0
    std_dict[i] = 0.0
    
#calculate average and standard deviation of all network parameters including number of wires and number of junctions
for i in avgs_dict.keys():
    for j in params_dict.keys():
        avgs_dict[i] = avgs_dict[i] + params_dict[j][i]
    avgs_dict[i] = avgs_dict[i]/len(params_dict.keys())
    for j in params_dict.keys():
        std_dict[i] = (params_dict[j][i]-avgs_dict[i])**2
    std_dict[i] = np.sqrt(std_dict[i]/len(params_dict.keys()))

print avgs_dict
print std_dict