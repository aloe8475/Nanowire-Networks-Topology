# Comparing Networks + Tasks:

# Nanowire Networks:

Two Tasks: Memory Capacity + Nonlinear Transformation

Test each of 300 Nanowire Networks and plot performance below:

- [-2,2] Volts input for MC task, Shortest Path Voltage input for Nonlinear Transformation Task

![ASN Networks logMC vs NLT](C:\Users\61424\Documents\GitHub\CODE\Data\Figures\Functional Connectivity\Network Level Comparisons\ASN Networks\ASN Networks logMC vs NLT.png)

## Network-Level Comparisons

What Network-Level Properties increase performance on each task - if any?

### Small Worldness

![Memory Capacity vs Nonlinear Transformation vs Small Worldness -All Networks](C:\Users\61424\Documents\GitHub\CODE\Data\Figures\Functional Connectivity\Network Level Comparisons\ASN Networks\Memory Capacity vs Nonlinear Transformation vs Small Worldness -All Networks.png)

### Avg Degree

![Memory Capacity vs Nonlinear Transformation vs Avg Degree - All Networks](C:\Users\61424\Documents\GitHub\CODE\Data\Figures\Functional Connectivity\Network Level Comparisons\ASN Networks\Memory Capacity vs Nonlinear Transformation vs Avg Degree - All Networks.png)

### Modularity

![Memory Capacity vs Nonlinear Transformation vs Modularity - All Networks](C:\Users\61424\Documents\GitHub\CODE\Data\Figures\Functional Connectivity\Network Level Comparisons\ASN Networks\Memory Capacity vs Nonlinear Transformation vs Modularity - All Networks.png)

---

### PCA

In order to try and separate out which network-level properties have an influence on task performance, we performed PCA analysis on the log MC vs NLT scatterplot, with num components = 2:

![ASN Networks PCA Component Correlations](C:\Users\61424\Documents\GitHub\CODE\Data\Figures\Functional Connectivity\Network Level Comparisons\ASN Networks\PCA\ASN Networks PCA Component Correlations.png)

Correlation Strengths (absolute value):

- r < 0.3 None or very weak 

- 0.3 < r <0.5 Weak 

- 0.5 < r < 0.7 Moderate

- r > 0.7 Strong

### T-Test

To back up the PCA analysis, we split the nanowire networks into the 50 highest and 50 lowest performing networks, for each task (each task has a unique set of 50 networks for each category)

We then computed independent sample t-Tests to compare the 50 highest and 50 lowest performing networks, on each Network-Level measure.

We corrected for multiple comparisons using the Bonferroni method with alpha = 0.05

#### NLT Task

![NanowireNetworks Ttest 50max vs 50min BarPlots NLT](C:\Users\61424\Documents\GitHub\CODE\Data\Figures\Functional Connectivity\Network Level Comparisons\ASN Networks\NanowireNetworks Ttest 50max vs 50min BarPlots NLT.png)

Note: 1star = p < 0.05, 2stars = p < 0.01, 3stars = p < 0.005, 4stars = p < 0.001

#### MC Task

![NanowireNetworks Ttest 50max vs 50min BarPlots MC](C:\Users\61424\Documents\GitHub\CODE\Data\Figures\Functional Connectivity\Network Level Comparisons\ASN Networks\NanowireNetworks Ttest 50max vs 50min BarPlots MC.png)

Note: 1star = p < 0.05, 2stars = p < 0.01, 3stars = p < 0.005, 4stars = p < 0.001

- For both tasks
  - Networks with negative-low degree skewness perform higher than those with high degree skewness
  - Networks with lower small worldness (M ~= 0.61) tend to perform better than networks with higher small worldness

- For the NLT task,
  - Networks with lower modularity (M ~= 0.57) tend to perform better than networks with higher modularity 
  - Networks with higher Mean Degree (M ~= 30 ) tend to perform better than networks with lower Mean Degree (M~=12)
  - These differences were not significant for the MC task

#### Network-Generating Parameters

I also extended this analysis to 3 parameters used to generate the networks: 

- Centroid Dispersion - How close to the center of the virtual plane the wires were dispersed, on average (higher score = closer to center)
- Mean Wire Length - How long each wire is, on average
- Std Wire Length - The standard deviation of wire lengths

![NanowireNetworks Parameters Ttest 50max vs 50min BarPlots NLT](C:\Users\61424\Documents\GitHub\CODE\Data\Figures\Functional Connectivity\Network Level Comparisons\ASN Networks\NanowireNetworks Parameters Ttest 50max vs 50min BarPlots NLT.png)

![NanowireNetworks Parameters Ttest 50max vs 50min BarPlots MC](C:\Users\61424\Documents\GitHub\CODE\Data\Figures\Functional Connectivity\Network Level Comparisons\ASN Networks\NanowireNetworks Parameters Ttest 50max vs 50min BarPlots MC.png)

- For both tasks
  - It seems that networks with wildly varying standard deviation in wire length (M ~= 50) perform much worse than those with more constant wire lengths (M~=11)
  - Networks with higher average wire lengths (M~=250) tend to perform better than those with lower average wire lengths (M~=150)
- For the Memory Capacity, networks that have higher centroid dispersion rates (M~= 275) tend to perform better than those who have lower rates (M~=225), however this was not significant for the NLT task 

---

### Delta

$$
\begin{equation}
NLT_1=\frac{NLT - NLT_{min}} {range(NLT) }
\tag{1}
\end{equation}
$$

$$
\begin{equation}
MC_1=\frac{MC - MC_{min}} {range(MC) }
\tag{2}
\end{equation}
$$

$$
\begin{equation}
\Delta_i=NLT_1 - MC_1
\tag{3}
\end{equation}
$$

Below graph shows NLT vs MC, performance for 300 NWNs, where the color is Delta. 

Networks with positive delta (yellow) are more likely to perform well on only NLT task, whereas networks with negative delta (purple) are more likely to perform well only on MC. Those with delta = 0 perform equally well (or badly) on both tasks.

<img src="Z:\Documents\CODE\Data\Figures\Functional Connectivity\Standardized MC vs Standardized NLT vs NLT minus MC.png" alt="Standardized MC vs Standardized NLT vs NLT minus MC" style="zoom:33%;" />

#### Delta Correlations:

I correlated Delta values for each network with Graph Theory measures and Parameter measures for that network, to see if there's any relationship between the measures and Delta:

![DeltaCorrelations - NLT and MC](C:\Users\61424\Documents\GitHub\CODE\Data\Figures\Functional Connectivity\Network Level Comparisons\ASN Networks\DeltaCorrelations - NLT and MC.PNG)

![NLT minus MC standardized accuracy correlation with PC and MZ](C:\Users\61424\Documents\GitHub\CODE\Data\Figures\Functional Connectivity\Node Level Comparisons\ASN Networks\NLT minus MC standardized accuracy correlation with PC and MZ.png)![NLT minus 2VMC standardized accuracy correlation with PC and MZ](C:\Users\61424\Documents\GitHub\CODE\Data\Figures\Functional Connectivity\Node Level Comparisons\ASN Networks\NLT minus 2VMC standardized accuracy correlation with PC and MZ.png)

This is a joint (2D) Histogram of Delta (equation 3) correlations between all 300 networks with Standardized PC and MZ in a 100x100 PC|MZ space

- For each pixel (100x100), correlate all networks. If the correlation is positive, it means more networks with a high Delta (i.e. Higher performance on NLT) have MZ/PC in those regions, whereas if correlation is negative more networks with low Delta (i.e. Higher Performance on MC) have MZ/PC in those regions.

### Mean Accuracy:

For each network, calculate a mean value between the standardised NLT and the standardised MC:

$$
\begin{equation}
{x_i}=\frac{NLT_1 +MC_1}{2}
\tag{4}
\end{equation}
$$

#### Mean Accuracy Correlations:

![MeanAccuracyCorrelations - NLT and MC](C:\Users\61424\Documents\GitHub\CODE\Data\Figures\Functional Connectivity\Network Level Comparisons\ASN Networks\MeanAccuracyCorrelations - NLT and MC.PNG)

![NLT and MC Mean accuracy correlations with PC and MZ - Less2Volts](C:\Users\61424\Documents\GitHub\CODE\Data\Figures\Functional Connectivity\Node Level Comparisons\ASN Networks\NLT and MC Mean accuracy correlations with PC and MZ - Less2Volts.png)

This is a joint (2D) Histogram of Mean Accuracy (equation 4) correlations between all 300 networks with Standardized PC and MZ in a 100x100 PC|MZ space

- For each pixel (100x100), correlate all networks. If the correlation is positive, it means networks perform better on both MC and NLT tasks, whereas if the correlation is negative, the network performs worse on both tasks. 

---

### **Repeated MC with [-2,2] volts as the input for each network:**

<img src="C:\Users\61424\Documents\GitHub\CODE\Data\Figures\Functional Connectivity\Network Level Comparisons\ASN Networks\Standardized 2V MC vs Standardized NLT vs NLT minus MC.png" alt="Standardized 2V MC vs Standardized NLT vs NLT minus MC" style="zoom:33%;" />

#### Delta:

![NLT minus 2VMC standardized accuracy correlation with PC and MZ](C:\Users\61424\Documents\GitHub\CODE\Data\Figures\Functional Connectivity\Node Level Comparisons\ASN Networks\NLT minus 2VMC standardized accuracy correlation with PC and MZ.png)

#### Mean Accuracy:

![Mean NLT and 2VMC accuracy correlations with PC and MZ](C:\Users\61424\Documents\GitHub\CODE\Data\Figures\Functional Connectivity\Node Level Comparisons\ASN Networks\Mean NLT and 2VMC accuracy correlations with PC and MZ.png)

---

## Node-Level Comparison:

Next, we compare node-level properties, namely PC and MZ, for all tasks to see whether there are clear delineations between networks that perform well on either/both tasks and those that do not.

#### MC + PCoeff:

- PCoeff - Measures how ‘well-distributed’ the links of a node are among different modules. The participation coefficient is close to 1 if its links are uniformly distributed among all the modules, and 0 if all its links are within its own module.
- MZ - Measures how ‘well-connected’ a node is to other nodes in the module. High values of indicate high within-module degrees and vice versa.

Plot all MC and PCoeff values of top 50 and bottom 50:

- Top 50 Networks = Red

- Bottom 50 Networks = Blue

NLT:

![GuimeraAmeral 50MinMax ALL Networks - NLT](C:\Users\61424\Documents\GitHub\CODE\Data\Figures\Functional Connectivity\Node Level Comparisons\ASN Networks\50 Min 50 Max\Vi (Shortest Path) Voltage\GuimeraAmeral 50MinMax ALL Networks - NLT.png)

MC:

- Shortest Path Voltage:



<img src="C:\Users\61424\Documents\GitHub\CODE\Data\Figures\Functional Connectivity\Node Level Comparisons\ASN Networks\50 Min 50 Max\Vi (Shortest Path) Voltage\GuimeraAmeral 50MinMax ALL Networks - MC.png" alt="GuimeraAmeral 50MinMax ALL Networks - MC" style="zoom:20%;" />

- [-2,2] Voltage for MC:

<img src="Z:\Documents\CODE\Data\Figures\Functional Connectivity\GuimeraAmeral 50MinMax ALL Networks - 3V MC.png" alt="GuimeraAmeral 50MinMax ALL Networks - 3V MC" style="zoom:%;" />

Histogram visualizations of the above:

- Shortest Path Voltage:

<img src="C:\Users\61424\Documents\GitHub\CODE\Data\Figures\Functional Connectivity\Node Level Comparisons\ASN Networks\50 Min 50 Max\Vi (Shortest Path) Voltage\PCoeff MZ Histograms 50MinMax ALL Networks - NLT + MC.png" alt="PCoeff MZ Histograms 50MinMax ALL Networks - NLT + MC" style="zoom:25%;" />

- Voltage of [-2,2] for MC:

<img src="Z:\Documents\CODE\Data\Figures\Functional Connectivity\PCoeff MZ Histograms 50MinMax - NLT + 2V MC.png" alt="PCoeff MZ Histograms 50MinMax - NLT + 2V MC" style="zoom: 25%;" />

Mean PC and MZ:

- Shortest Path Voltage:

![GuimeraAmeral 50MinMax Networks Avg - NLT + MC](C:\Users\61424\Documents\GitHub\CODE\Data\Figures\Functional Connectivity\Node Level Comparisons\ASN Networks\50 Min 50 Max\Vi (Shortest Path) Voltage\GuimeraAmeral 50MinMax Networks Avg - NLT + MC.png)

- Voltage of [-2,2] for MC:

![GuimeraAmeral 50MinMax Networks Avg - NLT + 2V MC](C:\Users\61424\Documents\GitHub\CODE\Data\Figures\Functional Connectivity\Node Level Comparisons\ASN Networks\50 Min 50 Max\Vi (Shortest Path) Voltage\GuimeraAmeral 50MinMax Networks Avg - NLT + 2V MC.png)

---

Correlations:

Lastly, we plot the Top 50 networks for each task, minus the Bottom 50 networks on the 100x100 Standardized PC|MZ plane.

NLT:

![Top50 minus Bottom50 standardized accuracy correlation with PC and MZ - NLT](C:\Users\61424\Documents\GitHub\CODE\Data\Figures\Functional Connectivity\Node Level Comparisons\ASN Networks\50 Min 50 Max\Vi (Shortest Path) Voltage\Top50 minus Bottom50 standardized accuracy correlation with PC and MZ - NLT.png)

Above, a zero score represents no difference in that pixel between the Top50 and the Bottom50 networks. 

- A positive score means that the Bottom networks are more likely to be in that space, whereas a negative score means Top networks are more likely to be in that space. 

MC:

![Top50 minus Bottom50 standardized accuracy correlation with PC and MZ - MC](C:\Users\61424\Documents\GitHub\CODE\Data\Figures\Functional Connectivity\Node Level Comparisons\ASN Networks\50 Min 50 Max\Vi (Shortest Path) Voltage\Top50 minus Bottom50 standardized accuracy correlation with PC and MZ - MC.png)

---

## Subgraph Analysis:

#### NLT:

![Top 15 Networks Subgraphs T = 5sec - NLT](Z:\Documents\CODE\Data\Figures\Functional Connectivity\Top 15 Networks Subgraphs T = 5sec - NLT.png)

![Bottom 15 Networks Subgraphs T = 5sec - NLT](Z:\Documents\CODE\Data\Figures\Functional Connectivity\Bottom 15 Networks Subgraphs T = 5sec - NLT.png)

<img src="C:\Users\61424\Documents\GitHub\CODE\Data\Figures\Functional Connectivity\Node Level Comparisons\ASN Networks\50 Min 50 Max\Vi (Shortest Path) Voltage\Subgraphs\GuimeraAmeral 50MinMax Networks Avg Subgraphs at T=5sec - NLT.png" alt="GuimeraAmeral 50MinMax Networks Avg Subgraphs at T=5sec - NLT" style="zoom:33%;" />

![PCoeff MZ Histograms 50MinMax Subgraph at T=5sec - NLT](C:\Users\61424\Documents\GitHub\CODE\Data\Figures\Functional Connectivity\Node Level Comparisons\ASN Networks\50 Min 50 Max\Vi (Shortest Path) Voltage\Subgraphs\PCoeff MZ Histograms 50MinMax Subgraph at T=5sec - NLT.png)

#### MC:

![Top 15 Networks Subgraphs At Max Conductance Timestep - 3V MC](Z:\Documents\CODE\Data\Figures\Functional Connectivity\Top 15 Networks Subgraphs At Max Conductance Timestep - 3V MC.png)

- Need to figure out a way to take subgraph without conductance thresholding.
- Can we link recurrence with MC / Structure? - HOW DO WE DO THIS?
- Look at Betweenness Centrality

---

# Beta Sweep:

MC vs NLT:

<img src="C:\Users\61424\Documents\GitHub\CODE\Data\Figures\Functional Connectivity\Network Level Comparisons\BetaSweep\Standardized MC vs Standardized NLT vs NLT minus MC - BetaSweep.png" alt="Standardized MC vs Standardized NLT vs NLT minus MC - BetaSweep" style="zoom:33%;" />

### Delta

#### Delta Correlations:

![DeltaCorrelations - Betasweep NLT and MC](C:\Users\61424\Documents\GitHub\CODE\Data\Figures\Functional Connectivity\Network Level Comparisons\BetaSweep\DeltaCorrelations - Betasweep NLT and MC.PNG)

![NLT minus MC standardized accuracy correlation with PC and MZ - BetaSweep](C:\Users\61424\Documents\GitHub\CODE\Data\Figures\Functional Connectivity\Node Level Comparisons\BetaSweep\NLT minus MC standardized accuracy correlation with PC and MZ - BetaSweep.png)

### Mean Accuracy

#### Mean Accuracy Correlations:

![MeanAccuracyCorrelations - Betasweep NLT and MC](C:\Users\61424\Documents\GitHub\CODE\Data\Figures\Functional Connectivity\Network Level Comparisons\BetaSweep\MeanAccuracyCorrelations - Betasweep NLT and MC.PNG)

![Mean NLT and MC accuracy correlations with PC and MZ - BetaSweep](C:\Users\61424\Documents\GitHub\CODE\Data\Figures\Functional Connectivity\Node Level Comparisons\BetaSweep\Mean NLT and MC accuracy correlations with PC and MZ - BetaSweep.png)

## Node-Level Comparison:

#### NLT:

![GuimeraAmeral 50MinMax Beta ALL Networks - NLT](C:\Users\61424\Documents\GitHub\CODE\Data\Figures\Functional Connectivity\Node Level Comparisons\BetaSweep\NLT\All Networks\GuimeraAmeral 50MinMax Beta ALL Networks - NLT.png)

#### MC:

![GuimeraAmeral 50MinMax Beta ALL Networks - MC](C:\Users\61424\Documents\GitHub\CODE\Data\Figures\Functional Connectivity\Node Level Comparisons\BetaSweep\MC\All Networks\GuimeraAmeral 50MinMax Beta ALL Networks - MC.png)

#### NLT & MC:

![GuimeraAmeral 50MinMax Networks Beta Avg - NLT + MC](C:\Users\61424\Documents\GitHub\CODE\Data\Figures\Functional Connectivity\Node Level Comparisons\BetaSweep\GuimeraAmeral 50MinMax Networks Beta Avg - NLT + MC.png)

![PCoeff MZ Histograms 50MinMax Beta ALL Networks - NLT + MC](C:\Users\61424\Documents\GitHub\CODE\Data\Figures\Functional Connectivity\Node Level Comparisons\BetaSweep\PCoeff MZ Histograms 50MinMax Beta ALL Networks - NLT + MC.png)



