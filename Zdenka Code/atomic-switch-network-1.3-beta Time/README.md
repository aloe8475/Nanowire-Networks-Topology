# atomic-switch-networks

This code generates the structure of networks of silver nanowires and
simulates the dynamics of the atomic switchs resistance, currents and
voltages.


## Installation

### Network generation in Python (optional)
Two python scientific libraries are required: 
+ `numpy` and 
+ `scipy` 

### 3D visualizations in Python (optional)
+ `mayavi` (which runs in python 2)

## Quick Start: Simulations in Matlab
The main code has been tested with Matlab R2015a and R2015b.

Examples should run out of the box.

The basic usage to get started with simulation is described below

##### 1- Make sure you add the folder `asn` and its subdirectories to the matlab path. Change to `examples` dicrectory. Most of the useful demos are there.

##### 2- Set some general simulation parameters (eg, how long the simulation is, what to plot at the end, etc)

```Matlab
    %% Plot and analysis output flags:
    SimulationOptions.takingSnapshots = true;  
    SimulationOptions.onlyGraphics    = false; 

    %% Simulation general options:
    SimulationOptions.dt = 1e-3;   % (sec) integration time step size
    SimulationOptions.T  = 6e1;    % (sec) duration of simulation
    SimulationOptions.TimeVector = (SimulationOptions.dt:SimulationOptions.dt:SimulationOptions.T)';
    SimulationOptions.NumberOfIterations = length(SimulationOptions.TimeVector);  
 ```   
##### 3- Specify and the load a connectivity matrix (eg, the physical network)

```Matlab
    Connectivity.WhichMatrix    = 'nanoWires';
    Connectivity.filename = '2016-09-08-155153_asn_nw_00100_nj_00261_seed_042_avl_100.00_disp_10.00.mat';             
    Connectivity = getConnectivity(Connectivity);
```        
##### 4- Specify the type of `Component` and its dynamics, and initialize their defualt values

```Matlab
    Components.ComponentType       = 'atomicSwitch'; % 'atomicSwitch' \ 'memristor' \ 'resistor'
    Components = initializeComponents(Connectivity.NumberOfEdges,Components);
```

##### 5 - Specify the external stimulus

```Matlab
    Stimulus.BiasType = 'DC';
    Stimulus.Amplitude = 1.5;                             % (Volt)
    Stimulus = getStimulus(Stimulus, SimulationOptions);
```    
##### 6 - Define the contact points between which the external source and tester resistor will be conected 

```
    %% Simulation external source and recording options:
    SimulationOptions.ContactMode  = 'preSet';   
    SimulationOptions.ContactNodes = [1, 2]; % only really required for preSet, other modes will overwrite this
```

##### 7 - Generate the coefficients matrices derived from Kirchoff's laws

```Matlab
    %% Get equations: (KCL + KVL):
    Equations = getEquations(Connectivity,SimulationOptions.ContactNodes);
```

##### 8 - Run the simulation 

```Matlab
    [Output, SimulationOptions, snapshots] = simulateNetwork(Equations, Components, Stimulus, SimulationOptions);
```

##### 9 - Perform some analysis and plot the results

```Matlab
    plotResults(Output.networkResistance,Output.networkCurrent,Stimulus);
```
Almost all the options available, including more elaborated graphics and movie compilation, are found in [`examples/MonolithicDemo.m`](https://github.com/pausz/atomic-switch-network/blob/master/examples/MonolithicDemo.m)

## Contributing

1. Fork it!
2. Create your feature branch: `git checkout -b my-new-feature`
3. Commit your changes: `git commit -am 'Add some feature'`
4. Push to the branch: `git push origin my-new-feature`
5. Submit a pull request :D

## Credits

[Ido Marcus](https://github.com/ido-marcus)

[Paula Sanz-Leon](https://github.com/pausz/)

## License

TODO: Write license
