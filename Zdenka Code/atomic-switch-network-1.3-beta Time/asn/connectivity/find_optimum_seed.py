import scipy.io
import analyse_network
import numpy as np
import os
import operator

os.chdir('connectivity_data\\network_data-23-2-2017') #not sure if this will work on an os other than windows

params_dict = {} #empty dict where we will store network params
print 'there are ' + str(len(os.listdir('.'))) + ' networks to analyse'

desired_param = 'number_of_wires' #choose parameter to search

#loop through each file name
for i in os.listdir('.'):
    #load file and tie each parameter to be maximised to its generated seed
    wires_dict = scipy.io.loadmat(i)
    seed = wires_dict['this_seed']
    params_dict[int(seed)] = wires_dict[desired_param] #replace key here for a different network parameter

best_seed = max(params_dict.iteritems(), key = operator.itemgetter(1))[0] 
print str(best_seed) + ' is the seed which maximises ' + desired_param + ' to a value of ' + str(float(params_dict[best_seed]))