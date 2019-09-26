import numpy as np
import scipy.io as sio
import networkx as nx
import time

from tqdm import tqdm
from copy import deepcopy

class simulation_options__:
    def __init__(self, dt=1e-3, T=1e1, interfaceContactWires = [73,30]):
        self.dt = dt
        self.T = T
        self.TimeVector = np.arange(0, T, dt)
        self.NumOfIterations = self.TimeVector.size
        self.interfaceContactWires = interfaceContactWires
        self.contactWires = [interfaceContactWires[0]-1, interfaceContactWires[1]-1]

class connectivity__:
    def __init__(self, filename):
        fullpath = 'connectivity/connectivity_data/' + filename
        matfile = sio.loadmat(fullpath, squeeze_me=True, struct_as_record=False)
        for key in matfile.keys():
            if key[0:2] != '__':
                setattr(self, key, matfile[key])
        self.numOfJunctions = self.number_of_junctions
        self.numOfWires = self.number_of_wires

class junctionState__:
    def __init__(self, NumOfJunctions):
        self.type = 'Atomic_Switch'
        self.voltage = np.zeros(NumOfJunctions+1)
        self.resistance = np.zeros(NumOfJunctions+1)
        self.resistance[-1] = 0.1
        self.onResistance = np.ones(NumOfJunctions+1)*1e4
        self.onResistance[-1] = 0.1
        self.offResistance = np.ones(NumOfJunctions+1)*1e7
        self.offResistance[-1] = 0.1

        self.filamentState = np.zeros(NumOfJunctions+1)
        self.OnOrOff = np.full(NumOfJunctions+1, False, dtype=bool)
        self.setVoltage = np.ones(NumOfJunctions+1)*1e-2
        self.resetVoltage = np.ones(NumOfJunctions+1)*1e-3
        self.critialFlux = np.ones(NumOfJunctions+1)*1e-1
        self.maxFlux = np.ones(NumOfJunctions+1)*1.5e-1

    def updateResistance(self):
        self.resistance[-1] = 0.1
        self.OnOrOff[0:-1] = abs(self.filamentState[0:-1]) >= self.critialFlux[0:-1]
        self.resistance[0:-1] = self.offResistance[0:-1] + \
                                (self.onResistance[0:-1]-self.offResistance[0:-1])*self.OnOrOff[0:-1]

    def updateJunctionState(self, dt):
        self.filamentState = self.filamentState + \
                            (abs(self.voltage) > self.setVoltage) *\
                            (abs(self.voltage) - self.setVoltage) *\
                            np.sign(self.voltage) * dt
        self.filamentState = self.filamentState - \
                            (abs(self.voltage) < self.resetVoltage) *\
                            (self.resetVoltage - abs(self.voltage)) *\
                            np.sign(self.filamentState) * dt * 10    
        maxPosition = np.where(abs(self.filamentState) > self.maxFlux)[0]
        self.filamentState[maxPosition] = np.sign(self.filamentState[maxPosition]) * \
                                            self.maxFlux[maxPosition]

class stimulus__:
    def __init__(self, biasType='DC', onTime=1, offTime=2, TimeVector=np.arange(0, 1e1, 1e-3), 
                onAmp=1.1, offAmp=0.005):
        onIndex = np.where((TimeVector >= onTime) & (TimeVector <= offTime))
        self.signal = np.ones(TimeVector.size) * offAmp
        self.signal[onIndex] = onAmp

class equations__:
    def __init__(self, connectivity, contactWires):
        """
        Need to generate a (E+1) x (E+1) coefficient matrix.
        First V-1 equations based on Kirchoff's current law. 
        Sum of currents in junctions on the sum nanowire should be zero.
        """

        V = connectivity.numOfWires
        E = connectivity.numOfJunctions
        edgeList = connectivity.edge_list

        KCLmat = np.zeros((V,E+1))
        for this_wire in range(V):
            junction_ind, flow = np.where(edgeList == this_wire)
            for i in range(len(junction_ind)):
                KCLmat[this_wire, junction_ind[i]] = (flow[i] - 0.5)*2 

        KCLmat[contactWires[0], -1] = 1
        KCLmat[contactWires[1], -1] = -1

        """
        E-V+2 equations based on Kirchoff's voltage law.
        E-V+1 equations generated from cycle basis of the network itself.
        Last equation includes the external source.
        When have multiple electrodes, just include more equations.
        """

        KVLmat = np.zeros((E-V+2,E+1))
        G = nx.from_numpy_array(connectivity.adj_matrix)
        cycles = nx.cycle_basis(G)
        extCycle = nx.shortest_path(G, contactWires[1], contactWires[0])
        extCycle.append(V+1)
        cycles.append(extCycle)
        edgeList = np.append(edgeList, np.array([[contactWires[1], V+1]]), axis = 0)
        row_count = 0
        for this_cycle in cycles:
            temp_cycle = np.array(this_cycle)
            temp_cycle = np.append(temp_cycle, this_cycle[0])

            from_vertex = temp_cycle[0:-1]
            to_vertex = temp_cycle[1:]
            direc = ((from_vertex > to_vertex)-0.5)*2
            
            for i in range(direc.size):
                if not ((from_vertex[i] == contactWires[0]) and (to_vertex[i] == V+1)):
                    edge_pos = np.where((edgeList[:,0]==from_vertex[i])&(edgeList[:,1]==to_vertex[i]))[0]
                    if edge_pos.size != 0:
                        edge_index = edge_pos[0]
                    else: 
                        edge_index = np.where((edgeList[:,1]==from_vertex[i])&(edgeList[:,0]==to_vertex[i]))[0][0]
                KVLmat[row_count, edge_index] = direc[i]
            row_count += 1

        self.KCL = KCLmat[0:-1,:]
        self.KVL = KVLmat

    pass    

def simulateNetwork(simulationOptions, connectivity, junctionState, 
                    stimulus, equations, 
                    simpleOutput = False,
                    useSparse = False):

    niterations = simulationOptions.NumOfIterations
    E = connectivity.numOfJunctions
    V = connectivity.numOfWires
    testerVoltage = np.zeros(niterations)
    rhs = np.zeros((E+1,1))
    if useSparse:
        from scipy.sparse import csc_matrix
        from scipy.sparse.linalg import spsolve

    if simpleOutput: 
        for i in tqdm(range(niterations), desc='Running Simulation '):
            junctionState.updateResistance()
            Rmat = np.ones((V-1,1))*junctionState.resistance
            lhs = np.vstack((equations.KCL/Rmat, equations.KVL)) 
            rhs[-1,0] = stimulus.signal[i]
            if useSparse:
                LHS = csc_matrix(lhs)
                RHS = csc_matrix(rhs)
                junctionState.voltage = spsolve(LHS, RHS)
            else: 
                junctionState.voltage = np.linalg.solve(lhs, rhs[:,0])
                
            junctionState.updateJunctionState(simulationOptions.dt)
            testerVoltage[i] = junctionState.voltage[-1]

        outDict = dict(
            testerVoltage = testerVoltage,
            networkCurrent = testerVoltage/junctionState.resistance[-1],
            networkResistance = (stimulus.signal/testerVoltage-1)*junctionState.resistance[-1]
            )
        return outDict
    
    else:
        import dataStruct 
        Network = dataStruct.network__()
        Network.filamentState = np.zeros((niterations, E))
        Network.junctionVoltage = np.zeros((niterations, E))
        Network.junctionResistance = np.zeros((niterations, E))
        Network.junctionSwitch = np.zeros((niterations, E), dtype = bool)

        for i in tqdm(range(niterations), desc='Running Simulation '):
            junctionState.updateResistance()
            Rmat = np.ones((V-1,1))*junctionState.resistance
            lhs = np.vstack((equations.KCL/Rmat, equations.KVL)) 
            rhs[-1,0] = stimulus.signal[i]
            if useSparse:
                LHS = csc_matrix(lhs)
                RHS = csc_matrix(rhs)
                junctionState.voltage = spsolve(LHS, RHS)
            else: 
                junctionState.voltage = np.linalg.solve(lhs, rhs[:,0])

            junctionState.updateJunctionState(simulationOptions.dt)
            testerVoltage[i] = junctionState.voltage[-1]

            Network.filamentState[i,:] = junctionState.filamentState[0:-1]
            Network.junctionVoltage[i,:] = junctionState.voltage[0:-1]
            Network.junctionResistance[i,:] = junctionState.resistance[0:-1]
            Network.junctionSwitch[i,:] = junctionState.OnOrOff[0:-1]

        Network.numOfWires = V
        Network.numOfJunctions = E
        Network.adjMat = connectivity.adj_matrix
        Network.contactWires = simulationOptions.interfaceContactWires
        Network.criticalFlux = junctionState.critialFlux[0]
        Network.TimeVector = simulationOptions.TimeVector
        Network.stimulus = stimulus.signal
        Network.networkCurrent = testerVoltage/junctionState.resistance[-1]
        Network.networkResistance = (stimulus.signal/testerVoltage-1)*junctionState.resistance[-1]
        Network.junctionList = np.add(connectivity.edge_list, 1).T
        Network.KVL = equations.KCL[:,0:-1]
        Network.KVL = equations.KVL[:,0:-1]
        # Network.LHS = lhs
        # Network.RHS = rhs
        Network.connectivity = connectivity
    return Network

def simulateNetwork_Lite(simulationOptions, 
                        connectivity, junctionState, 
                        stimulus, equations):
    niterations = simulationOptions.NumOfIterations
    E = connectivity.numOfJunctions
    V = connectivity.numOfWires
    testerVoltage = np.zeros(niterations)
    rhs = np.zeros((E+1,1))

    import dataStruct 
    Network = dataStruct.network__()
    Network.filamentState = np.zeros((niterations, E))
    Network.junctionVoltage = np.zeros((niterations, E))
    Network.junctionResistance = np.zeros((niterations, E))
    Network.junctionSwitch = np.zeros((niterations, E), dtype = bool)

    for i in tqdm(range(niterations), desc='Running Simulation '):
        junctionState.updateResistance()
        Rmat = np.ones((V-1,1))*junctionState.resistance
        lhs = np.vstack((equations.KCL/Rmat, equations.KVL)) 
        rhs[-1,0] = stimulus.signal[i]
        junctionState.voltage = np.linalg.solve(lhs, rhs[:,0])

        junctionState.updateJunctionState(simulationOptions.dt)
        testerVoltage[i] = junctionState.voltage[-1]

        Network.filamentState[i,:] = junctionState.filamentState[0:-1]
        Network.junctionVoltage[i,:] = junctionState.voltage[0:-1]
        Network.junctionResistance[i,:] = junctionState.resistance[0:-1]
        Network.junctionSwitch[i,:] = junctionState.OnOrOff[0:-1]

    Network.numOfWires = V
    Network.numOfJunctions = E
    Network.adjMat = connectivity.adj_matrix
    Network.contactWires = simulationOptions.interfaceContactWires
    Network.criticalFlux = junctionState.critialFlux[0]
    Network.TimeVector = simulationOptions.TimeVector
    Network.stimulus = stimulus.signal
    Network.networkCurrent = testerVoltage/junctionState.resistance[-1]
    Network.networkResistance = (stimulus.signal/testerVoltage-1)*junctionState.resistance[-1]
    Network.junctionList = np.add(connectivity.edge_list, 1).T
    Network.KVL = equations.KCL[:,0:-1]
    Network.KVL = equations.KVL[:,0:-1]
    # Network.LHS = lhs
    # Network.RHS = rhs
    Network.connectivity = connectivity

    return Network
