import numpy as np
import matplotlib.pyplot as plt
import sys
from utils import *
from dataStruct import *
"""
All interface indexing, like interfaceContactWires, use the index convention of MatLab,
which starts from 1.
"""
SimulationOptions = simulation_options__(dt = 1e-3, T = 3,
                                        interfaceContactWires = [73, 30])

Connectivity = connectivity__(
    filename = '2016-09-08-155153_asn_nw_00100_nj_00261_seed_042_avl_100.00_disp_10.00.mat')

JunctionState = junctionState__(Connectivity.numOfJunctions)
Stimulus = stimulus__(biasType = 'DC', 
                    TimeVector = SimulationOptions.TimeVector, 
                    onTime = 0, offTime = 1,
                    onAmp = 1.1, offAmp = 0.005)

Equations = equations__(Connectivity, SimulationOptions.contactWires)
this_realization = simulateNetwork(SimulationOptions, Connectivity, JunctionState, 
                                    Stimulus, Equations, 
                                    simpleOutput=False,
                                    useSparse=False)

this_realization.allocateData()
this_realization.draw(time = 0.9, JunctionsToObserve =[4, 229])

plt.style.use('classic')
#plt.plot(1/this_realization.networkResistance)



# sub1 = subnetwork_KVL__(this_realization, CenterJunctions = [4])

# sub1.draw(time = 0.9, JunctionsToObserve = [4,229])

# myloop = sub1.Loops[1]

# myloop.draw(time = 0.9, PathHighlight = this_realization.shortestPaths[0])