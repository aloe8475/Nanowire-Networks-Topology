# MatLab User Guide for ASN simulation
  * First, add the simulation folder with all subfolders into MatLab path (not necessary for now, will be useful when more functionalities are enabled). 
  * All-in-one simulation runner is ready here: https://github.com/rzhu40/ASN_simulation/blob/master/MatLab/runSimulation.m.

### If one hopes to set up the simulation step by step, here is the guide: 
  
  * For starters, one should set simulation options like total simulation time and time step size:

    ```matlab
    SimulationOptions.dt = 1e-3;   % (sec)
    SimulationOptions.T  = 1e0;    % (sec) duration of simulation
    SimulationOptions.TimeVector = (SimulationOptions.dt:SimulationOptions.dt:SimulationOptions.T)';
    SimulationOptions.NumberOfIterations = length(SimulationOptions.TimeVector);  
    ```

 * Then pick electrodes by nanowire index. Do **NOTICE** that the order matters.
  
    ```matlab
    SimulationOptions.ContactMode     = 'preSet'; 
    SimulationOptions.electrodes      = [73,30,83,88];
    SimulationOptions.numOfElectrodes = length(SimulationOptions.electrodes);
    ```

 * Next, one would have to load connecitivity information of the network (100 nanowire is an example, it's saved in this folder):
 
    ```matlab
    Connectivity.WhichMatrix = 'nanoWires';  
    Connectivity.filename = '2016-09-08-155153_asn_nw_00100_nj_00261_seed_042_avl_100.00_disp_10.00.mat';
    Connectivity = getConnectivity(Connectivity);
    ```
 
 * Then initialize network components:
 
    ```matlab
    Components.ComponentType       = 'atomicSwitch'; % 'atomicSwitch' \ 'memristor' \ 'resistor'
    Components = initializeComponents(Connectivity.NumberOfEdges,Components);
    ```
* The most important part is to initialize electrode conditions one-by-one in the order one set up before:
    * All electrodes have a so called *'signal'* (signal will be all zero if the electrode is a drain.).
    * Signals are stored in (Number of Electrodes x 1) cells named as **Signals**.
    * If the electrode is a source, specify input type, amplitude, etc.
      **Don't forget to store it in Signals**
    ```matlab
    Stimulus1.BiasType       = 'SinglePulse';           
    Stimulus1.OnTime         = 0.0; 
    Stimulus1.OffTime        = 1.0;
    Stimulus1.AmplitudeOn    = 1.4;
    Stimulus1.AmplitudeOff   = 0.005;
    Signals{1,1} = getStimulus(Stimulus1, SimulationOptions);
    ```
    * If the electrode is a drain, simply initialize by type drain, and store in Signals:
    ```matlab
    Stimulus2.BiasType       = 'Drain';           % 'DC' \ 'AC' \ 'DCandWait' \ 'Ramp'
    Signals{2,1} = getStimulus(Stimulus2, SimulationOptions);
    ```
* We save the space for the tedious work of setting up the last two electrodes. The last step is to run simulation:
    ```matlab
    [Output, SimulationOptions, snapshots] = simulateNetwork(Connectivity, Components, Signals, SimulationOptions);
    ```
* Output should be the struct one hopes to see simulation results. snapshots are more of the usage of visualization.
