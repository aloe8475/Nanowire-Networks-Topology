## Laplacian Spectral Plot: 

(https://www.frontiersin.org/articles/10.3389/fncom.2013.00189/full#h8)

![image-20200831102809223](C:\Users\61424\AppData\Roaming\Typora\typora-user-images\image-20200831102809223.png)

Normalised Laplacian:

![image-20200924112552032](C:\Users\61424\AppData\Roaming\Typora\typora-user-images\image-20200924112552032.png)

Spectral Plots:

![image-20200924112612559](C:\Users\61424\AppData\Roaming\Typora\typora-user-images\image-20200924112612559.png)

![img](https://www.frontiersin.org/files/Articles/60429/fncom-07-00189-HTML/image_m/fncom-07-00189-g001.jpg)

Results:

Given that similar spectra reflect similar network structures, it has been noted that the distribution of Laplacian eigenvalues can be used for classification of networks ([Ipsen and Mikhailov, 2001](https://www.frontiersin.org/articles/10.3389/fncom.2013.00189/full#B33); [Vukadinović et al., 2002](https://www.frontiersin.org/articles/10.3389/fncom.2013.00189/full#B65); [Banerjee and Jost, 2008b](https://www.frontiersin.org/articles/10.3389/fncom.2013.00189/full#B7); [Cetinkaya et al., 2012](https://www.frontiersin.org/articles/10.3389/fncom.2013.00189/full#B17))

**Smallest Eigenvalues**:

- The smallest eigenvalues of the Laplacian spectrum -meaning the first few eigenvalues of the labeled spectrum with 0 = λ1 ≤… ≤ λ*n* ≤ 2 – include information on the community structure of a network
  - **smaller eigenvalues indicate longer diffusion times, revealing a large proportion of intramodule connections and a low number of intermodule connections**
- The smallest non-zero eigenvalue λ2 (in connected networks) thus provides information on the best possible cut of the network into two modules ([Chung, 1996](https://www.frontiersin.org/articles/10.3389/fncom.2013.00189/full#B19)). 

**Largest Eigenvalues:**

- The largest eigenvalue of the Laplacian spectrum provides information on the level of “bipartiteness” of (subparts of) a network ([Bauer and Jost, 2012](https://www.frontiersin.org/articles/10.3389/fncom.2013.00189/full#B11))
- A (sub)network is said to be (fully) bipartite if its nodes can be divided into two groups in such a way that nodes within the same group are not connected. As a result, bipartite networks lack cycles with an odd number of nodes. **A division of nodes on basis of the positive and negative values in the largest eigenvector indicates the most bipartite configuration of the network.**
- 

**Largest EigenGap:**

- The divisions of all eigenvectors up to *vi* can subsequently be combined to separate the nodes of the network into *i* communities, where **a possible number of communities for the optimal division might be suggested by the largest eigengap (λ*i* + 1 − λ*i*)** ([Shi and Malik, 2000](https://www.frontiersin.org/articles/10.3389/fncom.2013.00189/full#B49); [Cheng and Shen, 2010](https://www.frontiersin.org/articles/10.3389/fncom.2013.00189/full#B18)).
- 

- **Directed vs Undirected**: First eigenvalues were smaller in the directed spectra compared to the first eigenvalues in the undirected spectra, suggests a stronger community structure in the directed networks 

