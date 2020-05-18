#!/usr/bin/env python
# coding: utf-8

# # Functional Comparison of ASN with WS and C. Elegans

# This notebook generates 24 different ASN networks with 300 nodes, each of which has differing parameters (Average Path length, Wire Dispersion and Centroid Dispersion). 
# 
# We then find the average degree for each ASN network, and use that as 2k to generate corresponding grid-like, small-world and random Watts-Strogatz networks.
# 
# We also load a sample C. Elegans network for comparison.

# In[30]:


get_ipython().system('jupyter notebook --version')
get_ipython().system('python --version')
get_ipython().system('conda --version')
get_ipython().system('jupyter trust RunTasksDifferentNetworks_MultipleNetworks.ipynb')


# In[31]:

cd

# In[32]:


cd "Documents/edamame"


# In[33]:


from scipy.io import loadmat, savemat
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import datetime
import networkx as nx
from edamame import *
from tqdm import tqdm_notebook
import os
import edamame.core.wires as wires
from random import choice
import warnings
from IPython.core.debugger import set_trace
import nct
import bct

#warnings.filterwarnings('ignore')


# ## Functions:

# In[34]:


import pickle 
import _pickle as cPickle
import gzip
def compressed_pickle(obj, filename,protocol=-1):
    with gzip.open(filename, 'wb') as f:
        cPickle.dump(obj, f, protocol)


# In[35]:


def decompress_pickle(file):
    with gzip.open(file, 'rb') as f:
        loaded_object = cPickle.load(f)
        return loaded_object


# In[36]:


# this_seed2=700
def connected_component_subgraphs(G):
    for c in nx.connected_components(G):
        yield G.subgraph(c)


# In[37]:


#Select Largest Components
def select_largest_component_new(wires_dict):
    """
    Find and select largest connected component of the original graph G.
    Throws away unconnected components and updates all the keys in wires_dict 
    """
#     def connected_component_subgraphs(G):
#         for c in nx.connected_components(G):
#             yield G.subgraph(c)
    
    wires_dict['G'] = max(connected_component_subgraphs(wires_dict['G']), key=len)
#     set_trace()
    nw = len(wires_dict['G'].nodes())
    nj = len(wires_dict['G'].edges())   
    
    logging.info("The largest component has %5d nodes and %6d edges", nw, nj)

    # Replace values in the dictionary
    wires_dict['number_of_wires']     = nw
    wires_dict['number_of_junctions'] = nj
    wires_dict['xa'] = wires_dict['xa'][wires_dict['G'].nodes()] 
    wires_dict['ya'] = wires_dict['ya'][wires_dict['G'].nodes()] 
    wires_dict['xb'] = wires_dict['xb'][wires_dict['G'].nodes()] 
    wires_dict['yb'] = wires_dict['yb'][wires_dict['G'].nodes()]
    wires_dict['xc'] = wires_dict['xc'][wires_dict['G'].nodes()] 
    wires_dict['yc'] = wires_dict['yc'][wires_dict['G'].nodes()]
 
    # Keep old edge_list
    old_edge_list = [(ii, kk) for ii, kk in  zip(wires_dict['edge_list'][:, 0], wires_dict['edge_list'][:, 1])]
    # Remove old edge list
    wires_dict = wires.remove_key(wires_dict, 'edge_list') 
    # Save indices of intersections in the old graph
    ind_dict = {key:value for value,key in enumerate(old_edge_list)}
    new_edge_list = sorted([kk if kk[0] < kk[1] else (kk[1], kk[0]) for kk in wires_dict['G'].edges()], key=lambda x: x[0])
    # Find intersection between the two sets
    inter = set(ind_dict).intersection(new_edge_list)
    # Retrieve edge indices/positions from the old list
    edges_idx = [ind_dict[idx] for idx in inter]
       
    # These have length equal to number of junctions -- only get the ones we need
    wires_dict['xi'] = wires_dict['xi'][edges_idx] 
    wires_dict['yi'] = wires_dict['yi'][edges_idx] 
    
    # Get contiguous numbering of nodes
    # Build node mapping 
    node_mapping    = {key:value for value, key in enumerate(sorted(wires_dict['G'].nodes()))}
    # This  step also renames edges list
    wires_dict['G'] =  nx.relabel_nodes(wires_dict['G'] , node_mapping)

    # Swap node vertices if vertex 0 is larger than vertex 1, then sort by first element
    wires_dict['edge_list'] = np.asarray(sorted([kk if kk[0] < kk[1] else (kk[1], kk[0]) for kk in wires_dict['G'].edges()], key=lambda x: x[0]))
    
    # Save adjacency matrix of new graph
    wires_dict = wires.remove_key(wires_dict, 'adj_matrix') 
    wires_dict = wires.generate_adj_matrix(wires_dict)

    wire_distances = wires.cdist(np.array([wires_dict['xc'], wires_dict['yc']]).T, np.array([wires_dict['xc'], wires_dict['yc']]).T, metric='euclidean')    
    wires_dict['wire_distances'] = wire_distances

    return wires_dict 


# ## Generate 300nw ASN + C. Elegans:

# In[38]:


cElegans=loadmat("../CODE/Data/Organic Networks Connectomes/celegans277neurons.mat")

    #C:\Users\aloe8475\Documents\PhD\GitHub\CODE\Data\Organic Networks Connectomes\


# In[39]:


#set up dictionary for celegans:
Elegans={'adj_matrix':[],'G':[],'Accuracy':{'Linear Transformation':[],'Mackey Glass':{'Cheat Steps':[],'Accuracy Value':[]},'MNIST':[]},'Graph Theory':{'Small World':[],'Modularity':[],'CCoeff':[],'MZ':[],'PCoeff':[],'PL':[]}}


# In[40]:


elegansMat=cElegans['celegans277matrix']
elegansGraph = nx.from_numpy_array(elegansMat)
Elegans['adj_matrix']=elegansMat
Elegans['G']=elegansGraph


# In[41]:


cd "/import/silo2/aloe8475/Documents/CODE/Analysis/Functional Connectivity/Functional Tasks/"


# In[43]:


if ( os.path.isfile('networks_LinearTransformation.pkl')): #if we haven't saved the file
    print('Creating Parameters')
    #Set Paramaters for network generation: Centroid Dispersion, Wire Dispersion + Length of Wires
    params={'Centroid':np.arange(150,400,25),'Wire Dispersion':[1.0, 2.5, 5.0, 10.0, 25.0, 50.0, 75.0, 80.0, 90.0, 100.0],'Length':np.arange(100,350,25)}
    numNetworks=len(params['Centroid'])+len(params['Wire Dispersion'])+len(params['Length']) #all parameters)
    if 'ASN300' in locals():
        loaded=True
    else:
        loaded=False
        ASN300=[[None]*10 for i in range(numNetworks)]
        cluster1=[[None]*10 for i in range(len(params['Centroid']))] #change in Centroid Dispersion
        cluster2=[[None]*10 for i in range(len(params['Wire Dispersion']))] #change in Wire Dispersion
        cluster3=[[None]*10 for i in range(len(params['Length']))] #change in average wire length


    #loop through these:
    centroid2=params['Centroid']
    length2=params['Length']
    disp2=params['Wire Dispersion']
    for i in range(49):
        temp=params['Centroid']
        temp2=params['Length']
        temp3=params['Wire Dispersion']
        centroid2=np.concatenate((centroid2,temp))
        length2=np.concatenate((length2,temp2))
        disp2=np.concatenate((disp2,temp3))
    print('Parameters Saved')
else:
    print('Parameters Loaded')


# In[44]:


#This loops through 30 different parameters, and for each parameter generates networks until we find 10 networks with 280+ nodes,
# for a total of 300 networks of 280+ nodes with different parameters


#I create 30 sets of networks. For each set, i generate networks, changing 2 parameters, 
# until i've found 10 networks with >280 nodes. 
#This way I generate 300 networks with varying parameters.

if ( os.path.isfile('networks_LinearTransformation.pkl')): #if we haven't saved the file
    print('Networks not loaded, creating now')
    count1=0
    count2=0
    count3=0
    this_seed=1779#np.random.randint(10000)
#     seed2=range(1, 6001, 20)
    for i in range(numNetworks):
        #loop through different parameter sets:
        x=0
        for j in range(500): #for each parameter, create 100 different networks
            if i < len(params['Centroid']):
#                 seed2=np.random.randint(100000)
                #Change centroid dispersion + avg wire length, but keep dispersion constant
                if x < 10:
                    print('Parameter: ' + str(i+1) + ' , Network: ' + str(j+1))
                    print('Parameter 1: Centroid ' +str(params['Centroid'][i])+ ' , Parameter 2: Avg Length ' + str(length2[j]))
                    temp=wires.generate_wires_distribution(300,this_seed=np.random.randint(100000),wire_dispersion=10,centroid_dispersion=params['Centroid'][i],wire_av_length=length2[j])
                    temp=wires.detect_junctions(temp)
                    temp=wires.generate_adj_matrix(temp)
                    temp=wires.generate_graph(temp)
                    temp=select_largest_component_new(temp)
                    counter1=0
                    
                while temp['G'].number_of_nodes() <277 and counter1<=10: #if nodes are less than 277, try the same parameters again for 10 more times before moving on to the next parameter
                    print(str(counter1) +': Parameter 1: Centroid ' +str(params['Centroid'][i])+ ' , Parameter 2: Avg Length ' + str(length2[j]))
                    temp=wires.generate_wires_distribution(300,this_seed=np.random.randint(100000),wire_dispersion=10,centroid_dispersion=params['Centroid'][i],wire_av_length=length2[j])
                    temp=wires.detect_junctions(temp)
                    temp=wires.generate_adj_matrix(temp)
                    temp=wires.generate_graph(temp)
                    temp=select_largest_component_new(temp)
                    counter1=counter1+1
                    
                if temp['G'].number_of_nodes()>=277: #only networks with more than 277 nodes
                    if x < 10: #only store 10 networks of each type
                        ASN300[i][x]=temp         
                        cluster1[count1][x]=temp    
                    x = x+1
                else:
                    print('saved networks ' + str(x))

                #only increase count if it's the last j value in the loop
                if j == 499:
                    count1=count1+1

            elif i >= len(params['Centroid']) and i < len(params['Centroid'])+len(params['Wire Dispersion']):
                #Change centroid dispersion + wire dispersion, keep avg wire length constant
                if x < 10:
                    print('Parameter: ' + str(i+1) + ' , Network: ' + str(j+1))
                    print('Parameter 1: Wire Disp ' +str(params['Wire Dispersion'][i-len(params['Centroid'])])+ ' , Parameter 2: Centroid ' + str(centroid2[j]))
                    temp=wires.generate_wires_distribution(300,this_seed=np.random.randint(100000),centroid_dispersion=centroid2[j],wire_dispersion=params['Wire Dispersion'][i-len(params['Centroid'])],wire_av_length=110)
                    temp=wires.detect_junctions(temp)
                    temp=wires.generate_adj_matrix(temp)
                    temp=wires.generate_graph(temp)
                    temp=select_largest_component_new(temp)
                    
                    counter2=0
                while temp['G'].number_of_nodes() <277 and counter2<=10:
                    print(str(counter2) +': Parameter 1: Wire Disp ' +str(params['Wire Dispersion'][i-len(params['Centroid'])])+ ' , Parameter 2: Centroid ' + str(centroid2[j]))
                    temp=wires.generate_wires_distribution(300,this_seed=np.random.randint(100000),centroid_dispersion=centroid2[j],wire_dispersion=params['Wire Dispersion'][i-len(params['Centroid'])],wire_av_length=110)
                    temp=wires.detect_junctions(temp)
                    temp=wires.generate_adj_matrix(temp)
                    temp=wires.generate_graph(temp)
                    temp=select_largest_component_new(temp)
                    counter2=counter2+1
                    
                if temp['G'].number_of_nodes()>=277: #only networks with more than 250 nodes
                    if x < 10: #only store 10 networks of each type
                        ASN300[i][x]=temp         
                        cluster2[count2][x]=temp
                    x = x+1
                else:
                    print('saved networks ' + str(x))

                #only increase count if it's the last j value in the loop
                if j == 499:
                    count2=count2+1        
            else:
               #Change  wire dispersion +  avg wire length, keep centroid dispersion constant
                if x < 10:
                    print('Parameter: ' + str(i+1) + ' , Network: ' + str(j+1))
                    print('Parameter 1 Avg Length:  ' +str(params['Length'][i-(len(params['Centroid'])+len(params['Wire Dispersion']))])+ ' , Parameter 2: Wire Disp ' + str(disp2[j]))
                    temp=wires.generate_wires_distribution(300,this_seed=np.random.randint(100000),centroid_dispersion=300,wire_dispersion=disp2[j],wire_av_length=params['Length'][i-(len(params['Centroid'])+len(params['Wire Dispersion']))])
                    temp=wires.detect_junctions(temp)
                    temp=wires.generate_adj_matrix(temp)
                    temp=wires.generate_graph(temp)
                    temp=select_largest_component_new(temp)
                    
                    counter3=0
                while temp['G'].number_of_nodes() <277  and counter3<=10:
                    print(str(counter3) +': Parameter 1 Avg Length:  ' +str(params['Length'][i-(len(params['Centroid'])+len(params['Wire Dispersion']))])+ ' , Parameter 2: Wire Disp ' + str(disp2[j]))
                    temp=wires.generate_wires_distribution(300,this_seed=np.random.randint(100000),centroid_dispersion=300,wire_dispersion=disp2[j],wire_av_length=params['Length'][i-(len(params['Centroid'])+len(params['Wire Dispersion']))])
                    temp=wires.detect_junctions(temp)
                    temp=wires.generate_adj_matrix(temp)
                    temp=wires.generate_graph(temp)
                    temp=select_largest_component_new(temp)
                    counter3=counter3+1
                    
                if temp['G'].number_of_nodes()>=277: #only networks with more than 250 nodes
                    if x < 10: #only store 10 networks of each type
                        ASN300[i][x]=temp   
                        cluster3[count3][x]=temp
                    x = x+1
                else:
                    print('saved networks ' + str(x))
    #             ASN300[i][j]=wires.generate_wires_distribution(300,this_seed=seed2[j],centroid_dispersion=300,wire_dispersion=params['Wire Dispersion'][j],wire_av_length=params['Length'][i-(len(params['Centroid'])+len(params['Wire Dispersion']))],)
    #             count = 0;
    #             cluster3[count3][j]=wires.generate_wires_distribution(300,this_seed=seed2[j],centroid_dispersion=300,wire_dispersion=params['Wire Dispersion'][j],wire_av_length=params['Length'][i-(len(params['Centroid'])+len(params['Wire Dispersion']))],)

                #only increase count if it's the last j value in the loop
                if j == 499:
                    count3=count3+1
                    
    #Save networks so we don't have to run this every time
    name='networks_LinearTransformation.pkl'
    with open(name, 'wb') as f:
        pickle.dump([ASN300,cluster1,cluster2,cluster3], f)
        
else: #load pickle file + communicability matrix calculated in Linear Transformation section 
    name='networks_LinearTransformation.pkl'
    print('Loading Networks + Linear Transformation Results')
    file = open(name, 'rb')
    [ASN300,cluster1,cluster2,cluster3,time_index,nodesList] = pickle.load(file)
#     [ASN300,cluster1,cluster2,cluster3] = pickle.load(file)
    print('Loaded')


# In[20]:


#set up accuracy dictionary for ASN

numberOfNodeTests=6
numberOfCheatSteps=6
count1=0
count2=0
count3=0
for i in range(len(ASN300)):
    for j in range(len(ASN300[i])):
        ASN300[i][j].update({'Accuracy':{'Linear Transformation':[None]*numberOfNodeTests,'Mackey Glass':{'Cheat Steps':[],'Accuracy Value':[None]*numberOfCheatSteps},'MNIST':[]},'Graph Theory':{'Small World':[],'Modularity':[],'CCoeff':[],'MZ':[],'PCoeff':[],'PL':[]}})
    #set up accuracy dictionary for ASN Clusters
        if i < len(cluster1):
            cluster1[count1][j].update({'Accuracy':{'Linear Transformation':[None]*numberOfNodeTests,'Mackey Glass':{'Cheat Steps':[],'Accuracy Value':[None]*numberOfCheatSteps},'MNIST':[]},'Graph Theory':{'Small World':[],'Modularity':[],'CCoeff':[],'MZ':[],'PCoeff':[],'PL':[]}})
            if j == len(ASN300[i]):
                count1=count1+1
        elif i >= len(cluster1) and i < len(cluster1)+len(cluster2):
            cluster2[count2][j].update({'Accuracy':{'Linear Transformation':[None]*numberOfNodeTests,'Mackey Glass':{'Cheat Steps':[],'Accuracy Value':[None]*numberOfCheatSteps},'MNIST':[]},'Graph Theory':{'Small World':[],'Modularity':[],'CCoeff':[],'MZ':[],'PCoeff':[],'PL':[]}})
            if j == len(ASN300[i]):
                count2=count2+1
        else:
            cluster3[count3][j].update({'Accuracy':{'Linear Transformation':[None]*numberOfNodeTests,'Mackey Glass':{'Cheat Steps':[],'Accuracy Value':[None]*numberOfCheatSteps},'MNIST':[]},'Graph Theory':{'Small World':[],'Modularity':[],'CCoeff':[],'MZ':[],'PCoeff':[],'PL':[]}})
            if j == len(ASN300[i]):
                count3=count3+1


# In[21]:


# # nwEdges100=np.array(np.sum(ASN100['adj_matrix'])/2)

# #Detect junctions, create adj matrix and graph, and store number of junctions
nwEdges300=[[None]*10 for i in range(len(ASN300))]
# count1=0
# count2=0
# count3=0
for i in range(len(ASN300)):
    for j in range(len(ASN300[i])):
#         print('Network ' + str(i) +', Iteration ' +str(j))
#         ASN300[i][j]=wires.detect_junctions(ASN300[i][j])
#         ASN300[i][j]=wires.generate_adj_matrix(ASN300[i][j])
#         ASN300[i][j]=wires.generate_graph(ASN300[i][j])
#         ASN300[i][j]=select_largest_component_new(ASN300[i][j])
        nwEdges300[i][j]=np.array(np.sum(ASN300[i][j]['adj_matrix'])/2)
#         if i < len(cluster1):
#             print('Cluster 1')
#             cluster1[count1][j]=ASN300[i][j]
#             if j == len(ASN300[i])-1:
#                 count1=count1+1
#         elif i >= len(cluster1) and i < (len(cluster2) + len(cluster1)):
#             print('Cluster 2')
#             cluster2[count2][j]= ASN300[i][j]
#             if j == len(ASN300[i])-1:
#                 count2=count2+1
#         elif i >= (len(cluster2) + len(cluster1)) and i < (len(cluster1)+len(cluster2)+len(cluster3)):
#             print('Cluster 3')
#             cluster3[count3][j]= ASN300[i][j]
#             if j == len(ASN300[i])-1:
#                 count3=count3+1


# In[ ]:

#export adj matrices to calculate small worldness in matlab:
adj_mats={"AdjMat":[[None]*len(ASN300[0]) for i in range(len(ASN300))]}
for i in range(len(ASN300)):
    for j in range(len(ASN300[i])):
        adj_mats['AdjMat'][i][j]=(ASN300[i][j]['adj_matrix'])
savemat('300nwASN_multipleNetworks.mat',adj_mats)
#C:\Users\aloe8475\Documents\PhD\GitHub\CODE\Analysis\Functional Connectivity\Functional Tasks\

# ## Graph Theory Measures:

# ### ASN

# In[24]:


# if (not os.path.isfile('networks_LinearTransformation.pkl')): #if we haven't saved the file
for i in range(len(ASN300)):
    for j in range(len(ASN300[i])):
        ASN300[i][j].update({'Graph Theory':{'Small World':[],'Modularity':[],'CCoeff':[],'MZ':[],'PCoeff':[],'PL':[]}})
for i in range(len(cluster1)):
        for j in range(len(cluster1[i])):
            cluster1[i][j].update({'Graph Theory':{'Small World':[],'Modularity':[],'CCoeff':[],'MZ':[],'PCoeff':[],'PL':[]}})
for i in range(len(cluster2)):
        for j in range(len(cluster2[i])):
            cluster2[i][j].update({'Graph Theory':{'Small World':[],'Modularity':[],'CCoeff':[],'MZ':[],'PCoeff':[],'PL':[]}})
for i in range(len(cluster3)):
        for j in range(len(cluster3[i])):
            cluster3[i][j].update({'Graph Theory':{'Small World':[],'Modularity':[],'CCoeff':[],'MZ':[],'PCoeff':[],'PL':[]}})


# In[25]:

# Small Worldness: 
# ------------------------------------
# CALCULATED IN MATLAB: smallworldness.m 
# found in C:\Users\aloe8475\Documents\PhD\GitHub\CODE\Analysis\Functional Connectivity\Functional Tasks
# ------------------------------------
temp=loadmat(r'300nwASN_multipleNetworks_smallworld.mat')
smallworld=temp['smallworld']
del temp

# In[27]:


# if (not os.path.isfile('networks_LinearTransformation.pkl')): #if we haven't saved the file
# Modularity, PCoeff, Small Worldness & MZ:

ci = []
pcoeff= []
mz= []
clustering = []
count1=0
count2=0
count3=0
for i in tqdm(range(len(ASN300))):
    for j in range(len(ASN300[i])):
        ci,q=nct.community_louvain(ASN300[i][j]['adj_matrix'])
        pcoeff=bct.participation_coef(ASN300[i][j]['adj_matrix'],ci)
        mi=bct.module_degree_zscore(ASN300[i][j]['adj_matrix'],ci)
        clustering=nx.clustering(ASN300[i][j]['G'])
        ASN300[i][j]['Graph Theory']['PL']=dict(nx.all_pairs_shortest_path_length(ASN300[i][j]['G']))
        ASN300[i][j]['Graph Theory']['Modularity']=ci
        ASN300[i][j]['Graph Theory']['PCoeff']=pcoeff
        ASN300[i][j]['Graph Theory']['MZ']=mz
        ASN300[i][j]['Graph Theory']['Small World']=smallworld[i][j]
        ASN300[i][j]['Graph Theory']['CCoeff']=clustering
        ASN300[i][j]['Graph Theory']['Degree']=nx.degree(ASN300[i][j]['G'])

        if i < len(cluster1):
            cluster1[count1][j]['Graph Theory']['PL']=ASN300[i][j]['Graph Theory']['PL']
            cluster1[count1][j]['Graph Theory']['Modularity']= ASN300[i][j]['Graph Theory']['Modularity']
            cluster1[count1][j]['Graph Theory']['PCoeff']=ASN300[i][j]['Graph Theory']['PCoeff']
            cluster1[count1][j]['Graph Theory']['MZ']=ASN300[i][j]['Graph Theory']['MZ']
            cluster1[count1][j]['Graph Theory']['Small World']=ASN300[i][j]['Graph Theory']['Small World']
            cluster1[count1][j]['Graph Theory']['CCoeff']= ASN300[i][j]['Graph Theory']['CCoeff']
            cluster1[count1][j]['Graph Theory']['Degree']= ASN300[i][j]['Graph Theory']['Degree']
            if j == len(ASN300[i])-1:
                count1=count1+1
        elif i >= len(cluster1) and i < (len(cluster1) + len(cluster2)):
            cluster2[count2][j]['Graph Theory']['PL']=ASN300[i][j]['Graph Theory']['PL']
            cluster2[count2][j]['Graph Theory']['Modularity']= ASN300[i][j]['Graph Theory']['Modularity']
            cluster2[count2][j]['Graph Theory']['PCoeff']=ASN300[i][j]['Graph Theory']['PCoeff']
            cluster2[count2][j]['Graph Theory']['MZ']=ASN300[i][j]['Graph Theory']['MZ']
            cluster2[count2][j]['Graph Theory']['Small World']=ASN300[i][j]['Graph Theory']['Small World']
            cluster2[count2][j]['Graph Theory']['CCoeff']= ASN300[i][j]['Graph Theory']['CCoeff']
            cluster2[count2][j]['Graph Theory']['Degree']= ASN300[i][j]['Graph Theory']['Degree']
            if j == len(ASN300[i])-1:
                count2=count2+1
        else:
            cluster3[count3][j]['Graph Theory']['PL']=ASN300[i][j]['Graph Theory']['PL']
            cluster3[count3][j]['Graph Theory']['Modularity']= ASN300[i][j]['Graph Theory']['Modularity']
            cluster3[count3][j]['Graph Theory']['PCoeff']=ASN300[i][j]['Graph Theory']['PCoeff']
            cluster3[count3][j]['Graph Theory']['MZ']=ASN300[i][j]['Graph Theory']['MZ']
            cluster3[count3][j]['Graph Theory']['Small World']=ASN300[i][j]['Graph Theory']['Small World']
            cluster3[count3][j]['Graph Theory']['CCoeff']= ASN300[i][j]['Graph Theory']['CCoeff']
            cluster3[count3][j]['Graph Theory']['Degree']= ASN300[i][j]['Graph Theory']['Degree']
            if j == len(ASN300[i])-1:
                count3=count3+1


# In[28]:


def community_layout(g, partition):
    """
    Compute the layout for a modular graph.


    Arguments:
    ----------
    g -- networkx.Graph or networkx.DiGraph instance
        graph to plot

    partition -- dict mapping int node -> int community
        graph partitions


    Returns:
    --------
    pos -- dict mapping int node -> (float x, float y)
        node positions

    """

    pos_communities = _position_communities(g, partition, scale=3.)

    pos_nodes = _position_nodes(g, partition, scale=1.)

    # combine positions
    pos = dict()
    for node in g.nodes():
        pos[node] = pos_communities[node] + pos_nodes[node]

    return pos

def _position_communities(g, partition, **kwargs):

    # create a weighted graph, in which each node corresponds to a community,
    # and each edge weight to the number of edges between communities
    between_community_edges = _find_between_community_edges(g, partition)

    communities = set(partition.values())
    hypergraph = nx.DiGraph()
    hypergraph.add_nodes_from(communities)
    for (ci, cj), edges in between_community_edges.items():
        hypergraph.add_edge(ci, cj, weight=len(edges))

    # find layout for communities
    pos_communities = nx.spring_layout(hypergraph, **kwargs)

    # set node positions to position of community
    pos = dict()
    for node, community in partition.items():
        pos[node] = pos_communities[community]

    return pos

def _find_between_community_edges(g, partition):

    edges = dict()

    for (ni, nj) in g.edges():
        ci = partition[ni]
        cj = partition[nj]

        if ci != cj:
            try:
                edges[(ci, cj)] += [(ni, nj)]
            except KeyError:
                edges[(ci, cj)] = [(ni, nj)]

    return edges

def _position_nodes(g, partition, **kwargs):
    """
    Positions nodes within communities.
    """

    communities = dict()
    for node, community in partition.items():
        try:
            communities[community] += [node]
        except KeyError:
            communities[community] = [node]

    pos = dict()
    for ci, nodes in communities.items():
        subgraph = g.subgraph(nodes)
        pos_subgraph = nx.spring_layout(subgraph, **kwargs)
        pos.update(pos_subgraph)

    return pos

# ### C. Elegans

# In[29]:

#Small world calculated on C Elegans Matrix in smallworld.m in MATLAB
temp=loadmat(r'cElegans_smallworld.mat')
smallworld_elegans=temp['cElegansSW'][0]
del temp

# In[31]:

# Modularity, PCoeff, Small Worldness & MZ:
ci = []
pcoeff= []
mz= []

ci,q=nct.community_louvain(elegansMat)
pcoeff=bct.participation_coef(elegansMat,ci)
mz=bct.module_degree_zscore(elegansMat,ci)
Elegans['Graph Theory']['MZ']=mz
Elegans['Graph Theory']['PCoeff']=pcoeff
Elegans['Graph Theory']['Modularity']=ci
Elegans['Graph Theory']['Small World']=smallworld_elegans
Elegans['Graph Theory']['PL']=dict(nx.all_pairs_shortest_path_length(elegansGraph))
Elegans['Graph Theory']['CCoeff']=nx.clustering(elegansGraph)
Elegans['Graph Theory']['Degree']=nx.degree(elegansGraph)


# # Task 1: Linear  Transformation

# In[32]:


#Regression
def NOKEVregression(target,absV): 
    inputx=np.vstack((np.ones(len(target)),absV)).T
    a1=np.linalg.lstsq(inputx,target)
    return a1


# In[33]:


#Subgraph AdjMat

#Threshold by conductance - when tunnelling becomes appreciable (offResistance * 10)

def getOnGraph(network, this_TimeStamp = 0):
    edgeList = network.connectivity.edge_list
    adjMat = np.zeros((network.numOfWires, network.numOfWires))
#     set_trace()
    adjMat[edgeList[:,0], edgeList[:,1]] = (1/network.junctionResistance[this_TimeStamp,:])>1e-06#network.junctionSwitch[this_TimeStamp,:] #CHANGE THIS TO CONDUCTANCE THRESHOLD?
    adjMat[edgeList[:,1], edgeList[:,0]] = (1/network.junctionResistance[this_TimeStamp,:])>1e-06#network.junctionSwitch[this_TimeStamp,:] #CHANGE THIS TO CONDUCTANCE THRESHOLD?
    onGraph = nx.from_numpy_array(adjMat)
    onGraph=nx.DiGraph.to_undirected(onGraph)
    
    return onGraph


# In[34]:


def getSubGraphComm(network, this_TimeStamp = 0):
    onGraph = getOnGraph(network, this_TimeStamp)
    components = [i for i in nx.connected_components(onGraph)]
    giant_component = components[np.argmax([len(i) for i in nx.connected_components(onGraph)])]
    nodes = list(giant_component)
    commMat = np.zeros((network.numOfWires, network.numOfWires))
    subComm = nx.communicability(onGraph.subgraph(giant_component))
    for i in nodes:
        for j in nodes:
            commMat[i,j] = subComm[i][j]
    return commMat


# In[35]:

#Communicability + Current Matrices
def commCurr(sim):
    startTime=500
    timeSteps=50
    endTime=1500
    time_index=[startTime,endTime,timeSteps]#]len(sim.junctionResistance),timeSteps]
    currMat=[None]*len(range(startTime,endTime,timeSteps))#startTime,len(sim.junctionResistance),timeSteps))
    nodesListFull=[None]*len(range(startTime,endTime,timeSteps))#startTime,len(sim.junctionResistance),timeSteps))
    commu_Mat=[None]*len(range(startTime,endTime,timeSteps))#startTime,len(sim.junctionResistance),timeSteps))
    new_currGraph=[None]*len(range(startTime,endTime,timeSteps))#startTime,len(sim.junctionResistance),timeSteps))
    count = 0
    for i in tqdm(range(startTime,endTime,timeSteps)):#startTime,len(sim.junctionResistance),timeSteps)): #for each timestep
        currMat[count] = np.zeros((sim.numOfWires,sim.numOfWires))
        edgeList = sim.connectivity.edge_list
        currMat[count][edgeList[:,0], edgeList[:,1]] = sim.junctionVoltage[i,:]/sim.junctionResistance[i,:] #-1,:
        currMat[count] = currMat[count] + currMat[count].T
        currGraph = nx.from_numpy_array(currMat[count])
        subGraph = getOnGraph(sim, this_TimeStamp=i)
        commu_Mat[count]=getSubGraphComm(sim, this_TimeStamp=i)
        
        components = [j for j in nx.connected_components(subGraph)] #all connected nodes in subgraph

        max_ind = np.argmax([len(j) for j in nx.connected_components(subGraph)])
        currGraph = nx.subgraph(currGraph, components[max_ind])
        new_currGraph[count] = currGraph
        currMat[count] = np.array(nx.adjacency_matrix(new_currGraph[count]).todense())

#         commu = nx.communicability(new_currGraph[count])
#         commu_Edges[count]=commu
#         subSize = len(currGraph)
        nodesList=list(currGraph.nodes)
        nodesListFull[count]=nodesList

#         commu_Mat[count] = np.array([commu[k][j] for k in nodesList for j in nodesList]).reshape(subSize,subSize)
        count = count+1
        
    return nodesListFull,commu_Mat, currMat, new_currGraph, time_index


# ## ASN Networks:

# In[36]:

onAmp=[[] for i in range(len(ASN300))]
shortestPath=[[None]*10 for i in range(len(ASN300))]
for i in range(len(ASN300)):
    for j in range(len(ASN300[i])):
        temp=getFarthestPairing(ASN300[i][j]['adj_matrix'])
        shortestPath[i][j]=nx.shortest_path_length(ASN300[i][j]['G'],temp[0],temp[1])
        onAmp[i].append(shortestPath[i][j]/5)

# In[42]:

#Run Simulations
dt = 1e-2
f=0.5
Time=5
period=1/f
   #Instantiate Variables       
stimulus=[[[] for i in range(len(ASN300))],[[] for i in range(len(ASN300))]]
accSqu=[]
#    accTri=[]
#    accSaw=[]
#    accDbl=[]
maxSqu=[]
#    maxTri=[]
#    maxSaw=[]
#    maxDbl=[]
dt = 1e-2
f=0.5
Time=5
#figure out voltage for each network:

#Choose Electrode Pattern
for i in range(len(onAmp)): #for each parameter
    for j in range(len(onAmp[i])): #for each network
        stimulus[0][i].append((stimulus__(biasType='AC',onAmp=onAmp[i][j],TimeVector=np.arange(0,Time,dt),f=f)))
        stimulus[1][i].append((stimulus__(biasType='Drain',TimeVector=np.arange(0,Time,dt)))) #we don't want this drain to be active during training

#Initialise Output Variables
period=[]
TimeVector=[]
voltage=[]
conductance=[]
switches=[]


#Run Simulations
count1=0
count2=0
count3=0
for i in range(0,len(ASN300)): #for each parameter
    results_ASN=[None]*10 #for i in range(len(ASN300))]
    for j in range(len(ASN300[i])): #for each network:
        print('Parameter ' + str(i+1), ', Network ' + str(j+1))
        #Run Simulations
    #     results=[]
        # Connectivity=connectivity__('700nw_14533junctions.mat')
        stimulus2 = [item for item in stimulus] #go through each list in the list and find the ith item
    #     set_trace()
        results_ASN[j]=runSim(connectivity__(wires_dict=ASN300[i][j]),stimulus=stimulus2, contactMode='farthest', T = Time, dt = 0.001, onAmp = onAmp[i][j], biasType='AC',f=f,junctionMode='tunneling')
        #wires_dict=newNetworkTest[chosenNetwork])
        results_ASN[j].frequency=f
        results_ASN[j].dt=0.001
        period=1/f

    name= r'/import/silo2/aloe8475/Documents/CODE/Data/Functional Connectivity/simulations_LinearTransformation_' + str(i)
    print('Saving 10 Simulations for Parameter ' +str(i+1))
    compressed_pickle([results_ASN],name)
    
#     del nwSqu,ResultSqu,TimeVector

# In[ ]:


#Run Regressions:
if ( os.path.isfile(r'/import/silo2/aloe8475/Documents/CODE/Data/Functional Connectivity/simulations_LinearTransformation_29')):
    nodesList=[[None]*10 for i in range(len(ASN300))]
    count1=0
    count2=0
    count3=0
    for i in range(len(ASN300)): #for each parameter
        print('Loading Parameter ' + str(i+1) +' - All Networks')
        check_memory()
        [results_ASN]=decompress_pickle(r'/import/silo2/aloe8475/Documents/CODE/Data/Functional Connectivity/simulations_LinearTransformation_'+str(i))
#         del temp
        print('Loaded Parameter ' + str(i+1) +' - All Networks')
        check_memory()
        for j in range(len(ASN300[i])): #for each network:
            nwSqu =[None]*numberOfNodeTests
            print('Regressing Parameter ' + str(i+1), ', Network ' + str(j+1))
            TimeVector=results_ASN[j].TimeVector
            voltage=results_ASN[j].wireVoltage
            conductance=results_ASN[j].conductance
            switches=results_ASN[j].junctionSwitch

            target1= (onAmp[i][j] * (-np.sign(TimeVector % period - period/2)))
    #             target2[j] = (4*onAmp[i]/period * abs((TimeVector-period/4) % period - period/2) - onAmp[i])
    #             target3[j] = (onAmp[i]/period * (TimeVector % period))
    #             target4[j] = (onAmp[i]*np.sin(4*np.pi*(1/period)*TimeVector))

            if len(ASN300[i][j]['G']) >= 250:
                nodesList[i][j]=[50,100,150,200,250,len(ASN300[i][j]['G'])]#range(50, len(ws300[i][0])+1,50)
            elif len(ASN300[i][j]['G']) >= 200 and len(ASN300[i][j]['G']) < 250:
                nodesList[i][j]=[50,100,150,200,len(ASN300[i][j]['G'])]#range(50, len(ws300[i][0])+1,50)  
            elif len(ASN300[i][j]['G']) >= 150 and len(ASN300[i][j]['G']) < 200:
                nodesList[i][j]=[50,100,150,len(ASN300[i][j]['G'])]#range(50, len(ws300[i][0])+1,50)  
            elif len(ASN300[i][j]['G']) >= 100 and len(ASN300[i][j]['G']) < 150:
                nodesList[i][j]=[50,100,len(ASN300[i][j]['G'])]#range(50, len(ws300[i][0])+1,50)

            countK=0
            for k in nodesList[i][j]: #loop through sets of nodes for regression
                print('Running Regression: ' + str(k) + ' nodes')
                ResultSqu=[]
                ResultSqu = nonLinearTrans(results_ASN[j],'Square',k, repeats=50) #simulation, type of signal, number of nodes to sample from, number of linear regression repetitions (take avg)

            # OLD WAY: NOKEVregression(target1[j],nwOutputs[j].T)[0]
            # OLD WAY: outputx=np.vstack((np.ones(len(target1[j])),nwOutputs[j].T)).T
            # temp=np.dot(outputx,ResultSqu)
            # MSE=np.mean((target1[j]-temp)**2)
            # rnMSE=np.sqrt(np.sum((target1[j]-temp)**2)/np.sum((target1[j])**2))
                nwSqu[countK]=ResultSqu['accuracy']
                ASN300[i][j]['Accuracy']['Linear Transformation'][countK]=nwSqu[countK]
                if i < len(cluster1):
                    cluster1[count1][j]['Accuracy']['Linear Transformation'][countK]=nwSqu[countK]
                    if j == len(ASN300[i]): #only increase count on the last 'j' network
                        count1=count1+1
                elif i >= len(cluster1) and i < (len(cluster1) + len(cluster2)):
                    cluster2[count2][j]['Accuracy']['Linear Transformation'][countK]=nwSqu[countK]
                    if j == len(ASN300[i]): #only increase count on the last 'j' network
                        count2=count2+1
                else:
                    cluster3[count3][j]['Accuracy']['Linear Transformation'][countK]=nwSqu[countK]
                    if j == len(ASN300[i]): #only increase count on the last 'j' network            
                        count3=count3+1
                countK=countK+1
                    #Save networks so we don't have to run this every time


# In[ ]:


#Current and Communicability Matrices:
if ( os.path.isfile(r'C:\Users\61424\Documents\GitHub\CODE\Data\Functional Connectivity\simulations_LinearTransformation_29')):
    time_index=[[None]*10 for i in range(len(ASN300))]
    junctions=[[None]*10 for i in range(len(ASN300))]
    sources=[[None]*10 for i in range(len(ASN300))]
    drains=[[None]*10 for i in range(len(ASN300))]
    nodesList=[[None]*10 for i in range(len(ASN300))]
    count1=0
    count2=0
    count3=0
    for i in tqdm(range(len(ASN300))):
        check_memory()
        #load corresponding results_ASN file
        print('Loading Parameter ' + str(i+1) +' - All Networks')
        [results_ASN,nodesList[i][j]]=decompress_pickle(r'C:\Users\61424\Documents\GitHub\CODE\Data\Functional Connectivity\simulations_LinearTransformation_'+str(i))
        print('Loaded Parameter ' + str(i+1) +' - All Networks')
        check_memory()
        for j in range(len(results_ASN)):
            print('Calculating COMM Matrix for Parameter ' + str(i+1) +', Network ' + str(j+1))
            nodesListASN_LT,commuMatASN_LT,currMatASN_LT,currGraphASN_LT,time_index[i][j]=commCurr(results_ASN[j]) #calculate communicability every 500 time steps for each network
            junctions[i][j]=results_ASN[j].junctionSwitch #save junction switch
            sources[i][j]=results_ASN[j].sources[0]
            drains[i][j]=results_ASN[j].drains[0]
            if i < len(cluster1):
                cluster1[count1][j]['Graph Theory']['COMM Mat']=commuMatASN_LT
                cluster1[count1][j]['Graph Theory']['Nodes List']=nodesListASN_LT
                cluster1[count1][j]['Graph Theory']['Current Matrix']=currMatASN_LT
                cluster1[count1][j]['Graph Theory']['Subgraph']=currGraphASN_LT
                if j == len(results_ASN)-1:
                    count1=count1+1
            elif i >= len(cluster1) and i < (len(cluster1) + len(cluster2)):
                cluster2[count2][j]['Graph Theory']['COMM Mat']=commuMatASN_LT
                cluster2[count2][j]['Graph Theory']['Nodes List']=nodesListASN_LT
                cluster2[count2][j]['Graph Theory']['Current Matrix']=currMatASN_LT
                cluster2[count2][j]['Graph Theory']['Subgraph']=currGraphASN_LT
                if j == len(results_ASN)-1: 
                    count2=count2+1
            else:
                cluster3[count3][j]['Graph Theory']['COMM Mat']=commuMatASN_LT
                cluster3[count3][j]['Graph Theory']['Nodes List']=nodesListASN_LT
                cluster3[count3][j]['Graph Theory']['Current Matrix']=currMatASN_LT
                cluster3[count3][j]['Graph Theory']['Subgraph']=currGraphASN_LT
                if j == len(results_ASN)-1:                
                    count3=count3+1
        del results_ASN, commuMatASN_LT, nodesListASN_LT, currMatASN_LT, currGraphASN_LT


# ## C. Elegans

# In[49]:


#On Amp as a function of shortestPath
onAmp=[]
shortestPath=[]
temp=getFarthestPairing(Elegans['adj_matrix'])
shortestPath=nx.shortest_path_length(Elegans['G'],temp[0],temp[1])
onAmp=shortestPath/5


# In[50]:


#Instantiate Variables       
stimulus_E=[[],[]]
accSqu=[]
#    accTri=[]
#    accSaw=[]
#    accDbl=[]
maxSqu=[]
#    maxTri=[]
#    maxSaw=[]
#    maxDbl=[]
dt = 1e-2
f=0.5
Time=5
# onAmp=2#[1.5,2,4,4,4,6,6,6,6,6,6,6,4,4,4,4,6,4,4,4,4,3,2,1.5]

#Choose Electrode Pattern

stimulus_E[0].append((stimulus__(biasType='AC',onAmp=onAmp,TimeVector=np.arange(0,Time,dt),f=f)))
stimulus_E[1].append((stimulus__(biasType='Drain',TimeVector=np.arange(0,Time,dt)))) #we don't want this drain to be active during training

#Initialise Output Variables
period=[]
TimeVector=[]
voltage=[]
conductance=[]
switches=[]
results=[None]*len([elegansGraph])
# Voltage=[None]*len(ASN300)
# Switches=[None]*len(ASN300)


nwSqu = []*len([elegansGraph])

#Run Simulations
# for i in range(len(ASN300)): #for each network
print('Network C.Elegans')
#Run Simulations
#     results=[]
# Connectivity=connectivity__('700nw_14533junctions.mat')
stimulus2 = [item for item in stimulus_E] #go through each list in the list and find the ith item
#     set_trace()
results=runSim(connectivity__(graph=elegansGraph),stimulus=stimulus2, contactMode='farthest', T = Time, dt = 0.001, onAmp = onAmp, biasType='AC',f=f,junctionMode='tunneling')
#wires_dict=newNetworkTest[chosenNetwork])
results.frequency=f
results.dt=0.001
period=1/f

TimeVector=results.TimeVector
voltage=results.wireVoltage
conductance=results.conductance
switches=results.junctionSwitch

stepNodes=len(elegansGraph)-1 #first use all nodes
sizes2=len(elegansGraph)
nwOutputs = [None]* int(sizes2/stepNodes)

outputNodes2=[]

for k in range(stepNodes,sizes2+1,stepNodes): #stepping up nodes
    np.random.seed(69)

    outputNodes2.append(voltage[:,np.random.choice(len(voltage[0,:]),size=k,replace=False)]) #take all the times (:) for a random j nodes

    nwOutputs=outputNodes2 #Length of nwOutputs is k (list)

target1=[None]*len(nwOutputs)
#     target2=[None]*len(nwOutputs)
#     target3=[None]*len(nwOutputs)
#     target4=[None]*len(nwOutputs)
#     target5=[None]*len(nwOutputs)

for j in range(len(nwOutputs)):
    target1[j] = (onAmp * (-np.sign(TimeVector % period - period/2)))
#         target2[j] = (4*onAmp[i]/period * abs((TimeVector-period/4) % period - period/2) - onAmp[i])
#         target3[j] = (onAmp[i]/period * (TimeVector % period))
#         target4[j] = (onAmp[i]*np.sin(4*np.pi*(1/period)*TimeVector))


#        nwTri     = np.zeros(len(nwOutputs))
#        nwSaw     = np.zeros(len(nwOutputs))
#        nwDbl     = np.zeros(len(nwOutputs))
#     nwMG      = [[None]*networksLoad,[None]*len(nwOutputs[0])]

ResultSqu=[]
#        ResultTri=[]
#        ResultSaw=[]
accuracy=[]
#        accuracyTri=[]
#        accuracySaw=[]
#        accuracyDbl=[]
output=[]
mSqu=[]
#        mTri=[]
#        mSaw=[]
#        mDbl=[]

# RUN LINEAR TRANSFORMATION
nodesListElegans=[50,100,150,200,250,len(Elegans['G'])]
for j in nodesListElegans: #range(50, len(Elegans['G'])+1,50): #loop through sets of nodes 
    print('Running Regression: ' + str(j) + ' nodes')
    
    ResultSqu=nonLinearTrans(results,'Square',j, repeats=50) #NOKEVregression(target1[j],nwOutputs[j].T)[0]
    
    nwSqu.append(ResultSqu['accuracy'])#(1 - rnMSE)
    
    
#THRESHOLD FOR COMMUNICABILITY MEASURE
results=[results]
# # From Matlab:
# resistance=results.junctionResistance
# for i in range(len(resistance)):
#     if resistance[i]<10000:
#         resistance[i][resistance[i]<10000]=1
# else:
#     resistance[i]=0


# In[51]:


#Current and Communicability Matrices
commuMatElegans_LT=[None]*len(results)
currMatElegans_LT=[None]*len(results)
# for i in tqdm(range(len(results))):
nodesListElegans_LT,commuMatElegans_LT,currMatElegans_LT,new_currGraphElegans_LT,time_indexElegans=commCurr(results[0])

Elegans['Graph Theory']['COMM Mat']=commuMatElegans_LT
Elegans['Graph Theory']['Nodes List']=nodesListElegans_LT
Elegans['Graph Theory']['Current Matrix']=currMatElegans_LT
Elegans['Graph Theory']['Subgraph']=new_currGraphElegans_LT


# In[52]:


Elegans['Accuracy']['Linear Transformation']=nwSqu


# In[53]:


Elegans['Accuracy']['Linear Transformation']

# In[57]:


ASNaccuracy=[[[None]*len(ASN300[0]) for j in range(len(ASN300))] for i in range(len(nodesList[0][-1][0]))]
cluster1accuracy=[[[None]*len(cluster1[0]) for j in range(len(cluster1))] for i in range(len(nodesList[0][-1][0]))]
cluster2accuracy=[[[None]*len(cluster2[0]) for j in range(len(cluster2))] for i in range(len(nodesList[0][-1][0]))]
cluster3accuracy=[[[None]*len(cluster3[0]) for j in range(len(cluster3))] for i in range(len(nodesList[0][-1][0]))]

for i in range(len(ASN300)): #for each parameter
    temp=[]
    for j in range(len(ASN300[i])): # for each network
        for k in range(len(ASN300[i][j]['Accuracy']['Linear Transformation'])):
            temp=ASN300[i][j]['Accuracy']['Linear Transformation'][k]
            if k == 0:
                ASNaccuracy[k][i][j]=temp
            elif k == 1:
                ASNaccuracy[k][i][j]=temp
            elif k == 2:
                ASNaccuracy[k][i][j]=temp
            elif k == 3:
                ASNaccuracy[k][i][j]=temp
            elif k == 4:
                ASNaccuracy[k][i][j]=temp
            elif k == 5:
                ASNaccuracy[k][i][j]=temp

for i in range(len(ASNaccuracy)): 
    count1=0
    count2=0
    count3=0
    for j in range(len(ASNaccuracy[i])):
        if j < len(cluster1):
            cluster1accuracy[i][count1]=ASNaccuracy[i][j]
            count1=count1+1
        elif j >= len(cluster1) and j < (len(cluster1) + len(cluster2)):
            cluster2accuracy[i][count2]=ASNaccuracy[i][j]
            count2=count2+1
        else:
            cluster3accuracy[i][count3]=ASNaccuracy[i][j]
            count3=count3+1


#Save Matrices For interactivy in Notebook:
#Update pickle file:
# if ( os.path.isfile(r'C:/Users/61424/Documents/GitHub/CODE/Analysis/Functional Connectivity/Functional Tasks/networks_LinearTransformation.pkl')):
name='networks_LinearTransformation.pkl'
with open(name, 'wb') as f:
    pickle.dump([ASN300,cluster1,cluster2,cluster3,time_index,nodesList], f)