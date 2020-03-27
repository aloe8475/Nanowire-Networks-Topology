import numpy as np
import networkx as nx

from tqdm import tqdm
from copy import deepcopy
import draw

class network__:
    pass

    def allocateData(self):
        print('Allocating data to obejcts....')
        self.gridSize = [self.connectivity.length_x, self.connectivity.length_y]
        self.junctionCurrent = self.junctionVoltage/self.junctionResistance
        # self.networkConductance = 1/self.networkResistance        
        self.isOnCurrentPath()
        # self.getWireVoltage()

        print('.')
        self.Junctions = [junction__(self, i) for i in range(0, self.numOfJunctions+1)]
        print('.')
        self.Nanowires = [nanowire__(self, i) for i in range(0, self.numOfWires+1)]
        print('.')
        self.JunctionsInd = np.sort([self.Junctions[i].index for i in range(1,len(self.Junctions))]).astype(int)
        self.WiresInd = np.sort([self.Nanowires[i].index for i in range(1,len(self.Nanowires))]).astype(int)
        
        print('Done!')

    def isOnCurrentPath(self):
        self.onCurrentPath = np.zeros((self.TimeVector.size, self.numOfJunctions), dtype = bool)
        edgeList = self.connectivity.edge_list
        contactWires = [self.contactWires[0]-1, self.contactWires[1]-1]

        for this_time in range(self.TimeVector.size):
            onMat = deepcopy(self.connectivity.adj_matrix)
            for i in range(self.numOfJunctions):
                if self.junctionSwitch[this_time,i] == 0:
                    onMat[edgeList[i,0], edgeList[i,1]] = 0
                    onMat[edgeList[i,1], edgeList[i,0]] = 0
            onGraph = nx.from_numpy_array(onMat)
            component = nx.node_connected_component(onGraph, contactWires[0])
            if contactWires[1] in component:
                for i in range(self.numOfJunctions):
                    if (self.junctionSwitch[this_time,i])&(edgeList[i,0] in component)&(edgeList[i,1] in component):
                        self.onCurrentPath[this_time,i] = True

    def draw(self, **kwargs):
        if not hasattr(self, 'Junctions'):
            self.allocateData()
            
        Lx, Ly = self.gridSize

        if 'TimeStamp' in kwargs:
            this_TimeStamp = kwargs['TimeStamp']
        elif 'time' in kwargs:
            if kwargs['time'] in self.TimeVector:
                this_TimeStamp = np.where(self.TimeVector == kwargs['time'])[0][0]
            elif (kwargs['time'] < min(self.TimeVector)) or (kwargs['time'] > max(self.TimeVector)):
                print('Input time exceeds simulation period.')
                this_TimeStamp = np.argmin(abs(self.TimeVector - kwargs['time']))
            else:
                this_TimeStamp = np.argmin(abs(self.TimeVector - kwargs['time']))
        else:
            this_TimeStamp = 0

        if 'JunctionsToObserve' in kwargs:
            JunctionsToObserve = kwargs['JunctionsToObserve']
        else:
            JunctionsToObserve = []

        if 'PathHighlight' in kwargs:
            PathHighlight = kwargs['PathHighlight']
        else:
            PathHighlight = []
        
        fig = draw.draw(self, PathHighlight = PathHighlight, 
                        JunctionsToObserve = JunctionsToObserve, 
                        TimeStamp = this_TimeStamp)
        return fig

class junction__:
    def __init__(self, network, index):
        self.index = index
        self.position = np.array([network.connectivity.xi[index-1], network.connectivity.yi[index-1]])
        self.filamentState = network.filamentState[:,index-1]
        self.switch = network.junctionSwitch[:,index-1]
        self.resistance = network.junctionResistance[:,index-1]
        self.voltage = network.junctionVoltage[:,index-1]
        self.current = network.junctionCurrent[:,index-1]
        self.contactWires = network.junctionList[:,index-1]
        self.onCurrentPath = network.onCurrentPath[:,index-1]
        self.conductance = 1/self.resistance
        
        # self.Node1Voltage = network.nodeVoltage[:,self.ContactNodes[0]-1]
        # self.Node2Voltage = network.nodeVoltage[:,self.ContactNodes[1]-1]

        if index == 0:
            for field in self.__dict__.keys():
                setattr(self, field, 0)

class nanowire__:
    def __init__(self, network, index):
        self.index = index
        self.centerPosition = np.array([network.connectivity.xc[index-1], network.connectivity.yc[index-1]])
        self.wireEnds =  np.array([network.connectivity.xa[index-1], network.connectivity.ya[index-1],
                                    network.connectivity.xb[index-1], network.connectivity.yb[index-1]])
        self.voltage = network.wireVoltage[:,index-1]
        self.adj = network.adjMat[:,index-1]
        self.contactWires = np.add(np.where(self.adj == 1), 1)[0]
        self.contactJunctions = np.add(np.where(network.junctionList[1,:] == index),1)
        self.contactJunctions = np.append(self.contactJunctions, np.add(np.where(network.junctionList[0,:] == index),1))
        self.onPathChecker = np.sum([network.Junctions[i].onCurrentPath for i in self.contactJunctions], axis = 0)

        """
        This part determines the current on this nanowire
        """
        if self.index != 0:
            xa, ya, xb, yb = self.wireEnds
            xc, yc = self.centerPosition
            
            xi, yi = network.connectivity.xi, network.connectivity.yi

            self.sectionCenterX = np.array([])
            self.sectionCenterY = np.array([])
            self.sectionCurrentX = np.array([])
            self.sectionCurrentY = np.array([])

            self.wireAngle = np.arctan2(yb-ya, xb-xa) % np.pi
            self.wireAngle = np.where(self.wireAngle <= np.pi/2, self.wireAngle, self.wireAngle - np.pi)

            if xa != xb:
                sortedInd = np.argsort([network.Junctions[i].position[0] for i in self.contactJunctions])
            else:
                sortedInd = np.argsort([network.Junctions[i].position[1] for i in self.contactJunctions])

            sortedContactJunctions = self.contactJunctions[sortedInd]
            sortedContactWires = self.contactWires[sortedInd]

            direction = ((sortedContactWires < index) - 0.5) * 2

            self.wireCurrents = np.zeros((network.TimeVector.size, self.contactJunctions.size-1))
            for i in range(network.TimeVector.size):
                self.wireCurrents[i,:] = np.cumsum(network.junctionCurrent[i,np.add(sortedContactJunctions[0:-1], -1)]*direction[0:-1])

            # self.sectionCenterX = np.append(self.sectionCenterX, np.mean([xi[np.add(sortedContactJunctions[0:-1], -1)], xi[np.add(sortedContactJunctions[1:], -1)]], axis = 0))
            # self.sectionCenterY = np.append(self.sectionCenterY, np.mean([yi[np.add(sortedContactJunctions[0:-1], -1)], yi[np.add(sortedContactJunctions[1:], -1)]], axis = 0))
            self.sectionCenterX = np.mean([xi[np.add(sortedContactJunctions[0:-1], -1)], xi[np.add(sortedContactJunctions[1:], -1)]], axis = 0)
            self.sectionCenterY = np.mean([yi[np.add(sortedContactJunctions[0:-1], -1)], yi[np.add(sortedContactJunctions[1:], -1)]], axis = 0)
            self.sectionCurrentX = np.cos(self.wireAngle) * self.wireCurrents
            self.sectionCurrentY = np.sin(self.wireAngle) * self.wireCurrents

            self.isElectrode = 0
            if (index in network.sources) or (index in network.drains):
                if index in network.sources:
                    self.isElectrode = 1
                else:
                    self.isElectrode = -1

                if xa != xb:
                    if xa < xb:
                        self.contactEnd = [xb, yb]
                    else :
                        self.contactEnd = [xa, ya]
                else:
                    if ya < yb:
                        self.contactEnd = [xb, yb]
                    else:
                        self.contactEnd = [xa, ya]

                self.totalCurrent = np.sum(network.junctionCurrent[:,np.add(sortedContactJunctions, -1)]*direction, axis = 1).reshape(network.TimeVector.size, 1)
                self.sectionCenterX = np.append(self.sectionCenterX, np.mean([xi[sortedContactJunctions[-1]-1], self.contactEnd[0]]))
                self.sectionCenterY = np.append(self.sectionCenterY, np.mean([yi[sortedContactJunctions[-1]-1], self.contactEnd[1]]))
                self.sectionCurrentX = np.append(self.sectionCurrentX, np.cos(self.wireAngle) * self.totalCurrent, axis = 1)
                self.sectionCurrentY = np.append(self.sectionCurrentY, np.sin(self.wireAngle) * self.totalCurrent, axis = 1)

        else:
            for field in self.__dict__.keys():
                setattr(self, field, 0)

class loop__:
    def __init__(self, network, mat1d, **kwargs):
        self.gridSize = network.gridSize
        self.TimeVector = network.TimeVector
        self.mat1d = mat1d

        if len(mat1d) == network.numOfJunctions:
            self.onLoopJunctionsInd = np.add(np.where(mat1d != 0)[0], 1).astype(int)
            self.Junctions = [network.Junctions[i] for i in self.onLoopJunctionsInd]
            self.onLoopJunctionsNum = len(self.Junctions)

            self.onLoopWiresInd = []
            for this_junction in self.Junctions:
                for i in this_junction.contactWires:
                    if not i in self.onLoopWiresInd:
                        self.onLoopWiresInd = np.append(self.onLoopWiresInd, i)
            
            self.onLoopWiresInd = self.onLoopWiresInd.astype(int)
            self.Nanowires = [network.Nanowires[i] for i in self.onLoopWiresInd]
            self.onLoopWiresNum = len(self.Nanowires)


        # when the input is a row of adjacency matrix.
        # make diagonal of adjacency matrix 1.
        elif len(mat1d) == network.numOfWires:
            self.onLoopWiresInd = np.add(np.where(mat1d != 0)[0], 1).astype(int)
            self.Nanowires = [network.Nanowires[i] for i in self.onLoopWiresInd]
            self.onLoopWiresNum = len(self.Nanowires)

            self.onLoopJunctionsInd = []
            for this_wire in self.Nanowires:
                for i in this_wire.contactJunctions:
                    if not i in self.onLoopJunctionsInd:
                        self.onLoopJunctionsInd = np.append(self.onLoopJunctionsInd, i)

            self.onLoopJunctionsInd = self.onLoopJunctionsInd.astype(int)
            self.Junctions = [network.Junctions[i] for i in self.onLoopJunctionsInd]
            self.onLoopJunctionsNum = len(self.Junctions)

        self.JunctionsInd = np.sort([self.Junctions[i].index for i in range(len(self.Junctions))]).astype(int)
        self.WiresInd = np.sort(self.onLoopWiresInd).astype(int)
    
    def draw(self, **kwargs):
        Lx, Ly = self.gridSize

        if 'TimeStamp' in kwargs:
            this_TimeStamp = kwargs['TimeStamp']
        elif 'time' in kwargs:
            if kwargs['time'] in self.TimeVector:
                this_TimeStamp = np.where(self.TimeVector == kwargs['time'])[0][0]
            elif (kwargs['time'] < min(self.TimeVector)) or (kwargs['time'] > max(self.TimeVector)):
                print('Input time exceeds simulation period.')
                this_TimeStamp = np.argmin(abs(self.TimeVector - kwargs['time']))
            else:
                this_TimeStamp = np.argmin(abs(self.TimeVector - kwargs['time']))
        else:
            this_TimeStamp = 0

        if 'JunctionsToObserve' in kwargs:
            JunctionsToObserve = kwargs['JunctionsToObserve']
        else:
            JunctionsToObserve = []

        if 'PathHighlight' in kwargs:
            PathHighlight = kwargs['PathHighlight']
        else:
            PathHighlight = []
        
        fig = draw.draw(self, PathHighlight = PathHighlight, 
                        JunctionsToObserve = JunctionsToObserve, 
                        TimeStamp = this_TimeStamp)
        return fig

class subnetwork_KVL__:

    def __init__(self, network, CenterJunctions = [], **kwargs):
        self.gridSize = network.gridSize
        self.TimeVector = network.TimeVector
        self.CenterJunctionsInd = CenterJunctions

        self.Loops = []
        self.Junctions = []
        self.Nanowires = []
        self.KVL_LoopsInd = np.array([])
        self.JunctionsInd = np.array([])
        self.WiresInd = np.array([])

        for this_junctionInd in self.CenterJunctionsInd:
            tempLoopsInd = np.where(network.KVL[:,this_junctionInd-1] != 0)[0]
            for i in tempLoopsInd:
                if i not in self.KVL_LoopsInd:
                    self.KVL_LoopsInd = np.append(self.KVL_LoopsInd, i)

        self.KVL_LoopsInd = np.sort(self.KVL_LoopsInd.astype(int))
        self.Loops.extend([loop__(network, network.KVL[i,:]) for i in self.KVL_LoopsInd])

        for this_loop in self.Loops:
            for this_wire in this_loop.Nanowires:
                if not this_wire in self.Nanowires:
                    self.Nanowires.append(this_wire)
                    self.WiresInd = np.append(self.WiresInd, this_wire.index)

            for this_junction in this_loop.Junctions:
                if not this_junction in self.Junctions:
                    self.Junctions.append(this_junction)
                    self.JunctionsInd = np.append(self.JunctionsInd, this_junction.index)
        
        self.WiresInd = np.sort(self.WiresInd).astype(int)
        self.JunctionsInd = np.sort(self.JunctionsInd).astype(int)

    def draw(self, **kwargs):
        Lx, Ly = self.gridSize

        if 'TimeStamp' in kwargs:
            this_TimeStamp = kwargs['TimeStamp']
        elif 'time' in kwargs:
            if kwargs['time'] in self.TimeVector:
                this_TimeStamp = np.where(self.TimeVector == kwargs['time'])[0][0]
            elif (kwargs['time'] < min(self.TimeVector)) or (kwargs['time'] > max(self.TimeVector)):
                print('Input time exceeds simulation period.')
                this_TimeStamp = np.argmin(abs(self.TimeVector - kwargs['time']))
            else:
                this_TimeStamp = np.argmin(abs(self.TimeVector - kwargs['time']))
        else:
            this_TimeStamp = 0

        if 'JunctionsToObserve' in kwargs:
            JunctionsToObserve = kwargs['JunctionsToObserve']
        else:
            JunctionsToObserve = []

        if 'PathHighlight' in kwargs:
            PathHighlight = kwargs['PathHighlight']
        else:
            PathHighlight = []
        
        fig = draw.draw(self, PathHighlight = PathHighlight, 
                        JunctionsToObserve = JunctionsToObserve, 
                        TimeStamp = this_TimeStamp)
        return fig