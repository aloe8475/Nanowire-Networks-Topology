# The Python User Guide
  * The python equivalent of Ido's MatLab code is here:
  https://github.com/rzhu40/ASN_simulation/tree/master/Python/Ido's%20method.
  * To use the improved method please go to:
  https://github.com/rzhu40/ASN_simulation/tree/master/Python/multiElec.
  * Again, all-in-one simulation runner is ready as: 
  https://github.com/rzhu40/ASN_simulation/blob/master/Python/multiElec/runMulti.py.
  
  * Step-by-step guide starts with importing all the libararies:
 
    ```python
    import numpy as np
    import matplotlib.pyplot as plt
    from utils import *
    from dataStruct import *
    ```
    
  * First set up simulation options:
    ```python
    SimulationOptions = simulation_options__(dt = 1e-3, T = 1,
                                        interfaceElectrodes = [73, 30])
    ```
  
  * Load connectivity information: 
    ```python
    Connectivity = connectivity__(
    filename = '2016-09-08-155153_asn_nw_00100_nj_00261_seed_042_avl_100.00_disp_10.00.mat')
    ```
    
  * Then set up electrode inputs(zeros all the time for drain signal), append to SimulationOptions.stimulus:
    ```python
    SimulationOptions.stimulus = []
    tempStimulus = stimulus__(biasType = 'DC', 
            TimeVector = SimulationOptions.TimeVector, 
            nanowire = SimulationOptions.interfaceElectrodes[0],
            onTime = 0, offTime = 1,
            onAmp = 1.4, offAmp = 0.005)
    SimulationOptions.stimulus.append(tempStimulus)

    tempStimulus = stimulus__(biasType = 'Drain', 
            TimeVector = SimulationOptions.TimeVector,
            nanowire = SimulationOptions.interfaceElectrodes[1])
    SimulationOptions.stimulus.append(tempStimulus)
    ```
   
  * Initialize junction states:
    ```python
    JunctionState = junctionState__(Connectivity.numOfJunctions)
    ```
  
  * Finally, one would be able to run simulation:
    ```python
    this_realization = simulateNetworkPlus(SimulationOptions, Connectivity, JunctionState)
    ```
  * The simulation result returns basic network dynamics in matrices, like junction voltage, resistance, filament state, network current, etc.
 #For further usage
  
  * One will get a lot more useful information about the network, and also junction/nanowire level data by running:
    ```python
    this_realization.allocateData()
    ```
  * To visualize the network, run:
    ```python
    this_realization.draw()
    ```
  * By default, it returns the network distribution at time = 0. One can also do something like this:
    ```python
    this_realization.draw(time = 0.9, 
                          PathHighlight = [1,2,3,4],
                          JunctionsToObserve = [6,7,8,9])                     
    ```
    
