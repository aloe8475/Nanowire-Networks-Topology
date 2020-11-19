# Beta Sweep

I wanted to compare our model's performance on WS networks with varying beta scores, to see how the SWP affects networks where the only difference is the re-wiring. I created 50 different parameters, and 50 networks per parameter:

- 10 WS beta values (0.1, 0.2 ... 1)

- 6 Avg Degrees (2, 4, 8, 16, 32)

- 50 networks per combination of parameters (5x10x50) for a total of 2500 networks

$V_i = \frac{\sigma{_{sd}}}{5}$ on a sinusoidal wave was used as the input voltage for all networks, and the nonlinear target was a square wave. 

## Average across 50 networks for each combination of parameters:

(e.g. Beta 0.1, avg deg of 2, take avg of 50 networks with those paramaters)

### Network Level Parameters :

Summary: Nonlinear Accuracy, Memory Capacity, Small Worldness, Modularity (q) and Avg Path Length

![BetaSweep Networks Parameters](C:\Users\61424\Documents\GitHub\CODE\Data\Figures\Functional Connectivity\Node Level Comparisons\BetaSweep\BetaSweep Networks Parameters.png)

---

**Nonlinear Transformation:**

![BetaSweep_NonlinearTransform_NetworkLevelMeasures_BetaValsGrouped](C:\Users\61424\Documents\GitHub\CODE\Data\Figures\Functional Connectivity\Network Level Comparisons\BetaSweep\NLT\Mean Across Networks Per Paramater\BetaSweep_NonlinearTransform_NetworkLevelMeasures_BetaValsGrouped.png)

There are no clear patterns defining beta's that are better or worse performing. There is also no clear pattern for small worldness / modularity scores influencing accuracy. The only clear pattern is that a high Avg Path length is associated with poorer performance.

![BetaSweep_NonlinearTransform_NetworkLevelMeasures](C:\Users\61424\Documents\GitHub\CODE\Data\Figures\Functional Connectivity\Network Level Comparisons\BetaSweep\NLT\Mean Across Networks Per Paramater\BetaSweep_NonlinearTransform_NetworkLevelMeasures.png)

Here it is clear that, on average, networks with Avg Degree of 16 perform better than networks with other Avg Degrees, and networks with Avg Deg of 2 perform much worse than others. 

**Memory capacity:**

![BetaSweep_MemoryCapacity_NetworkLevelMeasures_BetaValsGrouped](C:\Users\61424\Documents\GitHub\CODE\Data\Figures\Functional Connectivity\Network Level Comparisons\BetaSweep\MC\Mean Across Networks Per Paramater\BetaSweep_MemoryCapacity_NetworkLevelMeasures_BetaValsGrouped.png)

No clear patterns other than Avg Path Length, but this time the performance is not such a big difference.

![BetaSweep_MemoryCapacity_NetworkLevelMeasures](C:\Users\61424\Documents\GitHub\CODE\Data\Figures\Functional Connectivity\Network Level Comparisons\BetaSweep\MC\Mean Across Networks Per Paramater\BetaSweep_MemoryCapacity_NetworkLevelMeasures.png)

No clear patterns, except avg degree of 2 performs worse on average than other avg deg. 

---

### Node Level Parameters:

<u>**Unsorted**:</u> 

- Degree:

![Mean BetaSweep Networks Degree Curve Fitting Histograms](C:\Users\61424\Documents\GitHub\CODE\Data\Figures\Functional Connectivity\Node Level Comparisons\BetaSweep\Mean BetaSweep Networks Degree Curve Fitting Histograms.png)

- Betweenness Centrality:

![Mean BetaSweep Networks Betweenness Centrality Curve Fitting Histograms](C:\Users\61424\Documents\GitHub\CODE\Data\Figures\Functional Connectivity\Node Level Comparisons\BetaSweep\Mean BetaSweep Networks Betweenness Centrality Curve Fitting Histograms.png)

(Fixed Scale):

![Mean BetaSweep Networks Betweenness Centrality Curve Fitting Histograms Scale Fixed](C:\Users\61424\Documents\GitHub\CODE\Data\Figures\Functional Connectivity\Node Level Comparisons\BetaSweep\Mean BetaSweep Networks Betweenness Centrality Curve Fitting Histograms Scale Fixed.png)

betweenness centrality - 

- Within-Degree Module Z-Score:

![Mean BetaSweep Networks MZ Curve Fitting Histograms](C:\Users\61424\Documents\GitHub\CODE\Data\Figures\Functional Connectivity\Node Level Comparisons\BetaSweep\Mean BetaSweep Networks MZ Curve Fitting Histograms.png)

- Participation Coefficient:

![Mean BetaSweep Networks PCOEFF Curve Fitting Histograms](C:\Users\61424\Documents\GitHub\CODE\Data\Figures\Functional Connectivity\Node Level Comparisons\BetaSweep\Mean BetaSweep Networks PCOEFF Curve Fitting Histograms.png)

- Guimera Amaral Plane (MZ vs PCoeff)

![Mean BetaSweep_GuimeraAmaralPlane unsorted](C:\Users\61424\Documents\GitHub\CODE\Data\Figures\Functional Connectivity\Node Level Comparisons\BetaSweep\Mean BetaSweep_GuimeraAmaralPlane unsorted.png)



## All 2500 networks:

### Network Level Measures

Max 50 and Min 50 networks (Nonlinear transform)

![MaxMin 50 BetaSweep Networks](C:\Users\61424\Documents\GitHub\CODE\Data\Figures\Functional Connectivity\Network Level Comparisons\BetaSweep\NLT\All networks\MaxMin 50 BetaSweep Networks.png)

**Nonlinear Transformation:**

- Categorised by Avg Degree:

![BetaSweep_NonlinearTransform_NetworkLevelMeasures_AllNetworks](C:\Users\61424\Documents\GitHub\CODE\Data\Figures\Functional Connectivity\Network Level Comparisons\BetaSweep\NLT\All networks\BetaSweep_NonlinearTransform_NetworkLevelMeasures_AllNetworks.png)

- Categorised by Beta Value:

### ![BetaSweep_NonlinearTransform_NetworkLevelMeasures_BetaValsGrouped_AllNetworks](C:\Users\61424\Documents\GitHub\CODE\Data\Figures\Functional Connectivity\Network Level Comparisons\BetaSweep\NLT\All networks\BetaSweep_NonlinearTransform_NetworkLevelMeasures_BetaValsGrouped_AllNetworks.png)

**Memory Capacity**

- Categorised by Avg Degree:

![BetaSweep_MemoryCapacity_NetworkLevelMeasures_AllNetworks](C:\Users\61424\Documents\GitHub\CODE\Data\Figures\Functional Connectivity\Network Level Comparisons\BetaSweep\MC\All Networks\BetaSweep_MemoryCapacity_NetworkLevelMeasures_AllNetworks.png)

- Categorised by Beta Value:

![BetaSweep_MemoryCapacity_NetworkLevelMeasures_BetaValsGrouped_AllNetworks](C:\Users\61424\Documents\GitHub\CODE\Data\Figures\Functional Connectivity\Network Level Comparisons\BetaSweep\MC\All Networks\BetaSweep_MemoryCapacity_NetworkLevelMeasures_BetaValsGrouped_AllNetworks.png)

Next, I take the Max and Min networks across all Betas and Avg Degrees, and conduct T-Tests to see how they differ on Network-level Properties:

![BetaSweep Ttest 50max vs 50min BarPlots NLT](C:\Users\61424\Documents\GitHub\CODE\Data\Figures\Functional Connectivity\Network Level Comparisons\BetaSweep\NLT\All networks\BetaSweep Ttest 50max vs 50min BarPlots NLT.png)

![BetaSweep Ttest 50max vs 50min BarPlots MC](C:\Users\61424\Documents\GitHub\CODE\Data\Figures\Functional Connectivity\Network Level Comparisons\BetaSweep\MC\All Networks\BetaSweep Ttest 50max vs 50min BarPlots MC.png)

For both tasks we see similar trends for which networks tend to perform the best, and which perform the worst. It is clear that networks with Higher SW, higher Modularity, lower Beta Value and lower Avg Degree tend to perform better on both NLT and MC, and Degree Skewness has no impact. 

### Nodes Level Measures

It's clear that within a certain Beta Value, there is some difference between high performing networks, and low performing networks. That difference is partly explained by the Avg Degree and Avg Path Length, but there might be node-level explanations to this.

For this, I isolate only one particular Beta Value from the others, e.g. Beta = 0.1, to try and isolate what the differences are within the Beta Value between high performing and low perfoming networks:

Beta0.1:

![Beta01 Subset SW vs Modularity MC Scatter](C:\Users\61424\Documents\GitHub\CODE\Data\Figures\Functional Connectivity\Network Level Comparisons\BetaSweep\MC\All Networks\Beta Values\Beta01 Subset SW vs Modularity MC Scatter.png)

I then perform T-Tests to compare the highest and lowest performing networks (for each task) within Beta = 0.1 (for example) to see if there are any network-level differences:

![Beta01 Ttest 50max vs 50min BarPlots NLT](C:\Users\61424\Documents\GitHub\CODE\Data\Figures\Functional Connectivity\Network Level Comparisons\BetaSweep\NLT\All networks\Beta Values\Beta01 Ttest 50max vs 50min BarPlots NLT.png)

![ ](C:\Users\61424\Documents\GitHub\CODE\Data\Figures\Functional Connectivity\Network Level Comparisons\BetaSweep\MC\All Networks\Beta Values\Beta01 Ttest 50max vs 50min BarPlots MC.png)