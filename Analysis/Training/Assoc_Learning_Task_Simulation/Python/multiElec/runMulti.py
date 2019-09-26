import numpy as np
import matplotlib.pyplot as plt
import sys
from utils import *
from dataStruct import *
from draw_mpl import *
"""
All interface indexing, like interfaceContactWires, use the index convention of MatLab,
which starts from 1.
"""
SimulationOptions = simulation_options__(dt = 1e-3, T = 1,
                                        interfaceElectrodes = [73, 30, 88, 83])

# Connectivity = connectivity__(
#     filename = '2016-09-08-155153_asn_nw_00100_nj_00261_seed_042_avl_100.00_disp_10.00.mat')

Connectivity = connectivity__(
    filename = '2016-09-08-155044_asn_nw_00700_nj_14533_seed_042_avl_100.00_disp_10.00.mat')

SimulationOptions.stimulus = []
tempStimulus = stimulus__(biasType = 'DC', 
        TimeVector = SimulationOptions.TimeVector, 
        nanowire = SimulationOptions.interfaceElectrodes[0],
        onTime = 0, offTime = 1,
        onAmp = 1.1, offAmp = 0.005)
SimulationOptions.stimulus.append(tempStimulus)

tempStimulus = stimulus__(biasType = 'Drain', 
        TimeVector = SimulationOptions.TimeVector,
        nanowire = SimulationOptions.interfaceElectrodes[1])
SimulationOptions.stimulus.append(tempStimulus)

tempStimulus = stimulus__(biasType = 'DC', 
        TimeVector = SimulationOptions.TimeVector, 
        nanowire = SimulationOptions.interfaceElectrodes[2],
        onTime = 0, offTime = 1,
        onAmp = 1.4, offAmp = 0.005)
SimulationOptions.stimulus.append(tempStimulus)

tempStimulus = stimulus__(biasType = 'Drain', 
        TimeVector = SimulationOptions.TimeVector,
        nanowire = SimulationOptions.interfaceElectrodes[3])
SimulationOptions.stimulus.append(tempStimulus)

SimulationOptions.stimulus.append(tempStimulus)
JunctionState = junctionState__(Connectivity.numOfJunctions)

this_realization = simulateNetworkPlus(SimulationOptions, Connectivity, JunctionState)

draw_mpl(this_realization)