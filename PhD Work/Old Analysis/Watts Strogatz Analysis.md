## Watts-Strogatz Random and Lattice Networks:

## Network Level Analyses:	

### Accuracy

![50 Min Max Accuracy + Small Woldness - ASN vs WS ALL Networks](C:\Users\61424\Documents\GitHub\CODE\Data\Figures\Functional Connectivity\Network Level Comparisons\Nonlinear Transformation\50 Min Max Accuracy + Small Woldness - ASN vs WS ALL Networks.png)

### Degree

![50 Min Max Accuracy + Avg Degree - ASN vs WS ALL Networks](C:\Users\61424\Documents\GitHub\CODE\Data\Figures\Functional Connectivity\Network Level Comparisons\Nonlinear Transformation\50 Min Max Accuracy + Avg Degree - ASN vs WS ALL Networks.png)

### Modularity

![50 Min Max Accuracy + Modularity - ASN vs WS ALL Networks](C:\Users\61424\Documents\GitHub\CODE\Data\Figures\Functional Connectivity\Network Level Comparisons\Nonlinear Transformation\50 Min Max Accuracy + Modularity - ASN vs WS ALL Networks.png)



## Voltage Sweep:

I sorted the Nonlinear Transformation accuracies of random and grid networks from min to max accuracy (accuracies taken from using $V_i = \frac{\sigma{_{sd}}}{5}$).  Then I took the top 50 and bottom 50 performing networks, and performed a voltage sweep. 

### Random

![50Min50Max Accuracy Random Networks Voltage Sweep w Error Bar](C:\Users\61424\Documents\GitHub\CODE\Data\Figures\Functional Connectivity\Node Level Comparisons\50 Min 50 Max\50Min50Max Accuracy Random Networks Voltage Sweep w Error Bar.png)

### Grid

![50Min50Max Accuracy Grid Networks Voltage Sweep w Error Bar](C:\Users\61424\Documents\GitHub\CODE\Data\Figures\Functional Connectivity\Node Level Comparisons\50 Min 50 Max\50Min50Max Accuracy Grid Networks Voltage Sweep w Error Bar.png)

Degrees of freedom ==> Edges voltage increases due to filament state, edges on first current path will grow at the same time and similar rate, rest of network is low until first current path. Winner take all pathway.

- In random networks, there are multiple current paths (no winner take all) - meaning there are higher degrees of freedom.
- As soon as we introduce structure, we've immediately constrained the freedom of the network to perform. 

# Barabasi Scale-Free Networks

Lastly, I wanted to compare the model's performance on Scale-Free random networks, so we can have a complete representation of which structure is best for our model. 