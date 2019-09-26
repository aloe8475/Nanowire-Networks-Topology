import numpy as np
import networkx as nx
from copy import deepcopy
from itertools import islice

def findCurrent(network, numToFind = 1):
    edgeList = network.connectivity.edge_list
    contactWires = [network.contactWires[0]-1, network.contactWires[1]-1]
    numFound = 0
    Paths = []
    interfacePaths = []
    foundTime = []
    for this_time in range(network.TimeVector.size):
        onMat = deepcopy(network.connectivity.adj_matrix)
        for i in range(network.numOfJunctions):
            if network.junctionSwitch[this_time,i] == 0:
                onMat[edgeList[i,0], edgeList[i,1]] = 0
                onMat[edgeList[i,1], edgeList[i,0]] = 0
        onGraph = nx.from_numpy_array(onMat)
        tempPaths = [i for i in nx.all_simple_paths(onGraph, 72, 29)]
        pathFormed = len(tempPaths) > 0
        if pathFormed & (len(tempPaths) > numFound):
            for i in tempPaths:
                if i not in Paths:
                    Paths.append(i)
                    interfacePaths.append(np.add(i,1).tolist())
                    foundTime.append(network.TimeVector[this_time])
                    numFound+=1
        if numFound >= numToFind:
            break
    
    if numFound < numToFind:
        print(f'Unfortunately, only {numFound} current paths found in simulation time.')
    return interfacePaths, foundTime

def k_shortest_paths(network, k, weight=None, getLength=False):
    G = nx.from_numpy_array(network.adjMat)
    realSource = network.contactWires[0]-1
    realTarget = network.contactWires[1]-1
    realPathList = list((islice(nx.shortest_simple_paths(G, realSource, realTarget, weight=weight), k)))
    
    PathList = []
    LengthList = []
    for this_path in realPathList:
        PathList.append([i+1 for i in this_path])

    if getLength == False:
        return PathList
    else:
        for this_path in PathList:
            LengthList.append(get_path_length(G,this_path,weight))
        return PathList, LengthList

def get_path_length(G, path, weight=None):
    realPath = [i-1 for i in path]
    length = 0
    if weight == None:
        weight = 'weight'
    for i in range(len(realPath)-1):
        length += G[realPath[i]][realPath[i+1]][weight]
    return length

# def get_communicability(network, path)