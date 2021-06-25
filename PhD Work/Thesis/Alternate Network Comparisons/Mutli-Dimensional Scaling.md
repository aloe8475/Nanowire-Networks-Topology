

# Multi-Dimensional Scaling

For MDS, we need to create dissimilarity matrices, which compare each network to each other network.

To do so, I have tried 2 different methods: 

- GCD (based on https://www.nature.com/articles/srep04547)
- Portrait Divergence (based on https://doi.org/10.1007/s41109-019-0156-x)

## Graphlet-Based Network Comparison Distances (GCD)

from https://www.nature.com/articles/srep04547

**Graphlet**: A small induced subgraph of a large network that appears at any frequency and hence is independent of a null model

- An induced subgraph means that once you pick the nodes in the large network, you must pick all the edges between them to form the subgraph
- There are a total of 30 total graphlets that can be created (see below image) and applied to any graph:
  - For each graphlet, the authors use *automorphism orbits* to explain the relationship between nodes. For 4-node graphlets there are 15 total and 11 non-redundant *automorphism orbits*, and for 5-node graphlets there are 72 total and 56 non-redundant *automorphism orbits*.
  - The authors' findings suggest that using up to 4-node graphlets is the most efficient way to separate between networks

![figure1](https://media.springernature.com/lw685/springer-static/image/art%3A10.1038%2Fsrep04547/MediaObjects/41598_2014_Article_BFsrep04547_Fig1_HTML.jpg)

Results:

**Neuromorphic Neural Networks**:



**Beta Sweep WS Networks**:

E.g. Beta = 0.1 - Max and Min Networks:

![Beta 0.1 Max15Min15 Networks MultiDimensionalScaling - GCD](C:\Users\61424\Documents\GitHub\CODE\Data\Figures\Functional Connectivity\Multi Dimensional Scaling\GCD\Single Beta Values\Beta 0.1 Max15Min15 Networks MultiDimensionalScaling - GCD.png)

Beta = 0.5

![Beta 0.5 Max15Min15 Networks MultiDimensionalScaling - GCD](C:\Users\61424\Documents\GitHub\CODE\Data\Figures\Functional Connectivity\Multi Dimensional Scaling\GCD\Single Beta Values\Beta 0.5 Max15Min15 Networks MultiDimensionalScaling - GCD.png)

Beta = 1

![Beta 1 Max15Min15 Networks MultiDimensionalScaling - GCD](C:\Users\61424\Documents\GitHub\CODE\Data\Figures\Functional Connectivity\Multi Dimensional Scaling\GCD\Single Beta Values\Beta 1 Max15Min15 Networks MultiDimensionalScaling - GCD.png)



## Portrait Divergence

from https://doi.org/10.1007/s41109-019-0156-x

## Network Portraits

![figure1](https://media.springernature.com/lw685/springer-static/image/art%3A10.1007%2Fs41109-019-0156-x/MediaObjects/41109_2019_156_Fig1_HTML.png)

**Network Portraits** ([Bagrow et al. 2008](https://link.springer.com/article/10.1007/s41109-019-0156-x#ref-CR3)) - a graph invariant matrix representation of a network that is useful for visualization purposes but also capable of comparing pairs of networks.

- **Network Portrait** $B_{l,k}$ is an array with $k$ nodes at distance $l$ (shortest-path)

- $B_{0,k}$ stores the **number of nodes** in the graph

- $B_{1,k}$ stores the **Degree Distribution** = $P_k$ 
- $B_{2...N,k}$ stores the distribution of ***N*-nearest neighbors** for each node

- The **number of edges** $M$ is $\sum^N_{k=0}kB_{1,k}$

## Comparing networks by comparing portraits

**Network Portrait Divergences**, a principled information-theoretic measure for comparing networks, building graph-invariant distributions using the information contained within portraits.

The rows of *B* may be interpreted as probability distributions:
$$
P(k∣ℓ)=\frac{1}{N}B_{ℓ,k}
$$
This is the (empirical) **probability that a randomly chosen node will have *k* nodes at distance *ℓ***

1) Choose two nodes uniformly at random with replacement. The probability that they are connected is:
$$
\frac{∑_cn^2_c}{N_2},
$$
where $n_c$ is the number of nodes within connected component *c*, the sum $∑_cn^2_c$ runs over the number of connected components and $∑_cn_c=N$

2) **The probability the two nodes are at a distance *ℓ* from one another is:**
$$
\frac{paths_{length=l}}{paths}
$$
Number of paths of length *l* divided by number of total paths 

3) Lastly, the probability that one of the two nodes has *k*−1 other nodes at distance *ℓ* is given by:
$$
\frac{kB_{ℓ,k}}{∑_{k′}k′B_{ℓ,k′}}
$$
Combining all these = **a single distribution that encompasses the distances between nodes weighted by the “masses” or prevalence of other nodes at those distance**

- AKA - probability for choosing a pair of nodes at distance *ℓ* and for one of the two randomly chosen nodes to have *k* nodes at that distance *ℓ*

$$
P(k,ℓ)=\frac{∑_cn^2_c}{N_2}\frac{∑_{k′}k′B_{ℓ,k′}}{∑_cn^2_c}\frac{kB_{ℓ,k}}{∑_{k′}k′B_{ℓ,k′}}=\frac{kB_{ℓ,k}}{N^2}
$$

- this distribution is not normalized unless $G$ is connected



FINALLY: The **Network Portrait Divergence** *D*JS(*G*,*G*′) between two graphs *G* and *G*′ is the Jensen-Shannon divergence as follows:
$$
DJS(G,G′)≡\frac{1}{2}KL(P||M)+\frac{1}{2}KL(Q||M)
$$
Where $M=\frac{1}{2}(P+Q)M=\frac{1}{2}(P+Q)$  is the mixture distribution of *P* and *Q* 

P = first portrait, Q = second portrait 

## Conclusion

**The Network Portrait Divergence $0≤D_{JS}≤1$ provides a single value to quantify the dissimilarity of the two networks by means of their distance distributions, with smaller *D*JS for more similar networks and larger *D*JS for less similar networks.** 

- It compares networks based entirely on the structure of their respective topologies: the measure is independent of how the nodes in the network are indexed and, further, does not assume the networks are defined on the same set of nodes

## Results

![MultiDimensionalScaling - Portrait Divergence Nonlinear Transform](C:\Users\61424\Documents\GitHub\CODE\Data\Figures\Functional Connectivity\Multi Dimensional Scaling\Portrait Divergence\MultiDimensionalScaling - Portrait Divergence Nonlinear Transform.png)

