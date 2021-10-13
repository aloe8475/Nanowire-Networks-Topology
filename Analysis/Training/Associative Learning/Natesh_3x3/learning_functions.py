"Created on Thu Jul  1 13:16:32 2021 "
" @author: gnate "
# Define some noise functions
import numpy as np
import pandas as pd

def zero_bit_noise(N):
    noise_list = []
    noise_vec = np.zeros((1,N))
    noise_list.append(noise_vec)
    return noise_list

def one_bit_noise(N):
    noise_list = []
    for i in range(0,N):
        noise_vec = np.zeros((1,N))
        noise_vec[0,i] = 1
        noise_list.append(noise_vec)
    return noise_list

def xor_array(x,y):
    a = x.reshape((1,9))
    b = y.reshape((1,9))
    output_vec = np.zeros((1,a.shape[1]), dtype = 'uint')
    for i in range(0, a.shape[1]):
        output_vec[0, i] = a[0, i] ^ b[0, i]
    return output_vec     

def hamming_dist(x,y):
    dist = 0
    output_vec = xor_array(x,y)
    dist = sum(output_vec.T)
    return dist

## Main program
  
def generate_data():
    N = 9
    orig_inputs = ['111010111', '101010010', '101111101']
    orig_inputs_matrix = np.zeros((len(orig_inputs), 1, N), dtype = 'uint')
    for i in range(0,len(orig_inputs)):
        input_tuple = tuple(orig_inputs[i])
        for j in range(0, N):
            if input_tuple[j] == '1': orig_inputs_matrix[i,0,j] = 1
    # Generate 5x5 noisy data
    noise_vectors = []
    noise_vectors.append(zero_bit_noise(N))
    noise_vectors.append(one_bit_noise(N))
    #print('Noise vectors - ', noise_vectors)
    count_vec = np.zeros(2, dtype = 'uint')
    for i in range(0,2):
        count_vec[i] = np.math.factorial(N)/(np.math.factorial(i)*np.math.factorial(N-i))

    noisy_data_matrix = np.zeros((len(orig_inputs)*int(sum(count_vec)), 1, N), dtype = 'uint')
    row_count = 0
    for i in range(0,len(orig_inputs)):
        for j in range(0,2):
            bit_noise = np.array(noise_vectors[j], dtype = 'uint')
            for k in range(0, int(count_vec[j])):
                noisy_data_matrix[row_count] = xor_array(orig_inputs_matrix[i], bit_noise[k])
                row_count += 1

    final_data_matrix = np.zeros((len(orig_inputs)*int(sum(count_vec)), 1, N+1), dtype = 'uint')
    final_data_matrix[0:len(orig_inputs)*int(sum(count_vec)), 0, 0:N] = noisy_data_matrix.reshape(len(orig_inputs)*int(sum(count_vec)), N)

    for i in range(0,len(orig_inputs)):
        final_data_matrix[int(sum(count_vec))*i: int(sum(count_vec))*(i+1), 0, N] = i

    np.savetxt("3x3image_gen.csv", final_data_matrix.reshape(len(orig_inputs)*int(sum(count_vec)), N+1), delimiter=",")	
	
 
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
#        ax.view_init(azim=angle)"

# Load 3x3 dataset
def load_data():
    letters = pd.read_csv("./3x3image_gen.csv", sep=',', header = None)
    data = letters.values[:, 0:9]
    sample_num = 30      # Number of data samples
    
    # Target labels: 0 - 'z', 1 - 'v', 2 - 'n
    targets = letters.values[:, 9]
    inputs = np.zeros((sample_num,10))
    inputs[:,0:9] = data
    inputs[:,9] = np.ones(sample_num)   #Bias input
    
    # One-hot encoding of outputs
    onehot_outputs = create_target(targets) #1 0 0 = z, 0 1 0 = v, 0 0 1 = n
    
    # Standardize data - Uncomment to standardize data, not necessary
    # inputs = data - np.mean(data)
	# inputs = inputs/(np.std(data))"
    return inputs, onehot_outputs, sample_num, targets
    
    
    
# OTHER FUNCTIONS:
#Generating Electrode positions in Network - Written by Ruomin Zhu
def genGridNW(xa,xb,ya,yb,ex,ey):
  e = []
  for i in range(len(ex)):
      d = np.zeros(len(xa))
      for j in range(len(xa)):
          d[j]=dist((xa[j], ya[j]), (xb[j], yb[j]), (ex[i], ey[i]))
      e.append(np.argmin(d))
  return np.array(e)

def point_on_line(a, b, p):
  ap = p - a
  ab = b - a
  result = a + np.dot(ap, ab) / np.dot(ab, ab) * ab
  return result

def dist(p1,p2,p3):
  p1=np.array(p1) #xa ya
  p2=np.array(p2) #xb yb
  p3=np.array(p3) # ex ey (electrode placement)
  
  #determine whether closest point to electrode is on the line, or outside the line 
  t=point_on_line(p1,p2,p3)
  if t[1]<p1[1]:
      r = np.linalg.norm(p3-p1) #if point is outside left (xvalues)
  elif t[1] > p2[1]:
      r = np.linalg.norm(p3-p2) #if point is outside right (xvalues)
  else:
      r = np.linalg.norm(np.cross(p2-p1, p1-p3))/np.linalg.norm(p2-p1) #if point is inside xvalues
  return r
      
      
 #VISUALISE NETWORK STATE
def getWeightedGraph(edgeList,numWires,nwState):#, this_TimeStamp = 0):
#     edgeList = network['edge_list'],
    adjMat = np.zeros((numWires, numWires))
#     set_trace(),
    adjMat[edgeList[:,0], edgeList[:,1]] = nwState.voltage*nwState.conductance#network.junctionSwitch[this_TimeStamp,:] #CHANGE THIS TO CONDUCTANCE THRESHOLD?
    adjMat[edgeList[:,1], edgeList[:,0]] = nwState.voltage*nwState.conductance#network.junctionSwitch[this_TimeStamp,:] #CHANGE THIS TO CONDUCTANCE THRESHOLD?
    WeightedGraph = nx.from_numpy_array(adjMat)
    WeightedGraph=nx.DiGraph.to_undirected(WeightedGraph)
    
    return WeightedGraph

def draw_network_state(connectivity,nwState):
    adjMat=connectivity.adj_matrix
    graph=nx.from_numpy_array(adjMat)
    OGgraph=graph.copy()
    pos=nx.kamada_kawai_layout(OGgraph)
    numWires=graph.number_of_nodes()
    edgeList=np.array(list(graph.edges()))
    weightedSubGraph=getWeightedGraph(edgeList,numWires,nwState)
    minWeights=-1e-5 #currents
    maxWeights=1e-5 #currents
    
    #draw network
    #%matplotlib inline
    f,ax=plt.subplots(figsize=(10,6))
    G=weightedSubGraph
    edge_weights=nx.get_edge_attributes(G,'weight')
#     G.remove_edges_from((e for e, w in edge_weights.items() if w <1e-5)) ,
    edges=G.edges()
    weights=[G[u][v]['weight'] for u,v in edges]
    
    #draw OG graph,
    pos=nx.kamada_kawai_layout(OGgraph)
    h=nx.draw_networkx_nodes(OGgraph,pos=pos,node_color='grey',node_size=10,ax=ax)
    h.set_zorder(1),
    
    h2=nx.draw_networkx_edges(G,pos=pos,ax=ax,edge_color=weights,edge_cmap=plt.cm.bwr,edge_vmin=minWeights,edge_vmax=maxWeights)
#     if h2:,
#         h2.set_norm(clrs.SymLogNorm(10)),
#         h2.set_zorder(3),
    #             if j == 10 and i == 6:,
    #                 plt.colorbar(h2),
    nx.draw_networkx_nodes(G,pos=pos,nodelist=sources,node_color='#3f9b0b',node_size=120,node_shape ='*',ax=ax)
    nx.draw_networkx_nodes(G,pos=pos,nodelist=drain_pool,node_color='r',node_size=120,node_shape ='*',ax=ax)
    plt.colorbar(h2,label='current',ax=ax)
    plt.show()