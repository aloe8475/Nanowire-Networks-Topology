function [OutputDynamics, SimulationOptions, snapshots] = simulateNetworkAdrian(Equations, Components, Stimulus, SimulationOptions, Connectivity, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulates the network and finds the resistance between the two contacts
% as a function of time.
%
% ARGUMENTS: 
% Equations - Structure that contains the (abstract) matrix of coefficients
%             (as documented in getEquations) and the number of nodes in
%             the circuit.
% Components - Structure that contains the component properties. Every 
%              field is a (E+1)x1 vector. The extra component is the
%              tester resistor connected in series to the voltage and to 
%              the network.
% Stimulus - Structure that contains the details of the external stimulus
%            (time axis details, voltage signal).
% SimulationOptions - Structure that contains general simulation details that are indepedent of 
%           the other structures (eg, dt and simulation length);
% varargin - if not empty, contains an array of indidces in which a
%            snapshot of the resistances and voltages in the network is
%            requested. This indices are based on the length of the simulation.
% OUTPUT:
% OutputDynamics -- is a struct with the activity of the network
%                    .networkResistance - the resistance of the network (between the two 
%                     contacts) as a function of time.
%                    .networkCurrent - the overall current from contact (1) to contact (2) as a
%                     function of time.
% Simulationoptions -- same struct as input, with updated field names
% snapshots - a cell array of structs, holding the resistance and voltage 
%             values in the network, at the requested time-stamps.
        
% REQUIRES:
% updateComponentResistance
% updateComponentState
%
% USAGE:
%{
    Connectivity = getConnectivity(Connectivity);
    contact      = [1,2];
    Equations    = getEquations(Connectivity,contact);
    Components   = initializeComponents(Connectivity.NumberOfEdges,Components)
    Stimulus     = getStimulus(Stimulus);
    
    OutputDynamics = runSimulation(Equations, Components, Stimulus);
%}
%
% Authors:
% Ido Marcus
% Paula Sanz-Leon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    %Need to put adjacency matrix in 'Equations.AdjMat'

    %% Initialize:
    compPtr       = ComponentsPtr(Components);        % using this matlab-style pointer to pass the Components structure by reference
    niterations   = SimulationOptions.NumberOfIterations; 
    testerVoltage = zeros(niterations, 1);                      % memory allocation for the voltage on the tester resistor as function of time
    RHSZeros      = zeros(Equations.NumberOfEdges,1); % the first E entries in the RHS vector.
    avLambda      = zeros(niterations,1);
    lambda_vals   = zeros(niterations,Equations.NumberOfEdges);
    voltage_vals  = zeros(niterations,Equations.NumberOfEdges);
    OutputDynamics=[];
    %% Use sparse matrices:
%     Equations.KCLCoeff = sparse(Equations.KCLCoeff); %We don't calculate
%     these anymore
%     Equations.KVLCoeff = sparse(Equations.KVLCoeff);%We don't calculate
%     these anymore
    RHSZeros           = sparse(RHSZeros);
    
    %% If snapshots are requested, allocate memory for them:
    if ~isempty(varargin)
        snapshots           = cell(size(varargin{1}));
        snapshots_idx       = sort(varargin{1}); 
    else
        nsnapshots          = 10;
        snapshots           = cell(nsnapshots,1);
        snapshots_idx       = ceil(logspace(log10(1), log10(niterations), nsnapshots));
    end
    kk = 1; % Counter
    
    %% Solve equation systems for every time step and update:
    for ii = 1 : niterations
        % Show progress:
        progressBar(ii,niterations);     
        
        % Update resistance values:
        updateComponentResistance(compPtr); 
        
        idx = SimulationOptions.ContactNodes;

        %Here we express Kirchoffs Laws to solve for absolute voltages
        % fill extended matrix with electrodes
        %Change to ./

        Gmat = double(Equations.AdjMat) ./ Components.resistance(1:end,ones(Equations.NumberOfNodes,1));
        imat = Gmat; %conductance matrix
        imat = imat - diag(sum(imat,2)); %removing nanowires 'junctions' with themselves
        ifix = sparse(length(imat)+length(idx),1);
        %need to modify Stimulus.Signal(t,electrode) so it includes all
        %sources and drains
        for i=1:length(idx)
            ifix(Equations.NumberOfNodes+i) = Stimulus.Signal(ii,i); %voltage of the sources
            %All elements touching the electrode
            touching = nonzeros(Equations.AdjMat(idx(i),:).* (1:Equations.NumberOfNodes));
            imat(end+1,touching) = 1; %setting all wires touching electrodes to 1
            imat(idx(i),end+1) = 1; %add extra column for all the electrodes
%             imat(Electrodes.PosIndex(idx(i)),end+1) = 1; %add extra column for all the electrodes
        end

        %solve
        vsols=imat\ifix;

        % retrieve currents + voltages:
        [rows,cols,~] = find(Equations.AdjMat);
        Curr = Gmat;
        voltdif = Gmat;
        for i=1:length(rows)
            Curr(cols(i),rows(i)) = (vsols(cols(i))-vsols(rows(i)))...
                *Curr(cols(i), rows(i));
            voltdif(cols(i),rows(i))=(vsols(cols(i))-vsols(rows(i)));
        end

        simout.SumRule=sum(sum(Curr,2));
        k=1;
%         for i=1:Numel
%             if isequal(OpenFlag(i),1)
%                 simout.(strcat('I',Electrodes.Name{i}))=NaN;
%                 simout.(strcat('V',Electrodes.Name{i}))=NaN;        
%                 continue;
%             end
%             simout.(strcat('I',Electrodes.Name{i}))=vsols(length(Curr)+k);
%             simout.(strcat('V',Electrodes.Name{i}))=Electrodes.Value{i}(TimeInd);
%             k=k+1;
%         end
        simout.Currents=Curr;
        simout.Voltages=vsols(1:end-length(idx));
        simout.VoltDif=voltdif;
        compPtr.comp.voltage = voltdif;
        
        % Update element fields:
        updateComponentState(compPtr, Stimulus.dt);    % ZK: changed to allow retrieval of local values
        
        for j = 1:size(Connectivity.EdgeList,2)
            CurrentVec(j)    = Curr(Connectivity.EdgeList(1,j), Connectivity.EdgeList(2,j));
            lambdaVec(j)     = compPtr.comp.filamentState(Connectivity.EdgeList(1,j), Connectivity.EdgeList(2,j));   
            resistanceVec(j) = compPtr.comp.resistance(Connectivity.EdgeList(1,j), Connectivity.EdgeList(2,j));
            voltageVec(j)    = voltdif(Connectivity.EdgeList(1,j), Connectivity.EdgeList(2,j)); 
            onOrOffVec(j)    = compPtr.comp.OnOrOff(Connectivity.EdgeList(1,j), Connectivity.EdgeList(2,j));
        end
        
%         [lambda_vals(ii,:), voltage_vals(ii,:)] = updateComponentState(compPtr, Stimulus.dt);
        
        % Record tester voltage:
%         testerVoltage(ii) = compPtr.comp.voltage(end);
        
        %Record average value of lambda 
%         avLambda(ii) = mean(abs(compPtr.comp.filamentState));            
        
        % Record the activity of the whole network
        if find(snapshots_idx == ii) 
                frame.Timestamp  = SimulationOptions.TimeVector(ii);
                frame.Voltage    = voltageVec;
                frame.Currents   = CurrentVec;
                frame.Resistance = resistanceVec;
                frame.OnOrOff    = onOrOffVec;
                frame.filamentState = lambdaVec;
                frame.netV = Stimulus.Signal(ii);
%                 frame.netI = testerVoltage(ii) * compPtr.comp.resistance(end);
%                 frame.netC = 1/((Stimulus.Signal(ii) / testerVoltage(ii) - 1)/compPtr.comp.resistance(end));
                snapshots{kk} = frame;
                kk = kk + 1;
        end
        
            OutputDynamics.Currents{ii}=simout.Currents;
        
    end
    
    % Store some important fields
    SimulationOptions.SnapshotsIdx = snapshots_idx; % Save these to access the right time from .TimeVector.

    % Calculate network resistance and save:
%     OutputDynamics.testerVoltage     = testerVoltage;
%     OutputDynamics.networkCurrent    = testerVoltage .* compPtr.comp.resistance(end);
   %sum(OutputDynamics.Currents) is current flowing in out each node
   % take nth entry if nth is source / drain to get the source / drain
   % currents
   
    OutputDynamics.sumCurrents=sum(OutputDynamics.Currents{ii});
    OutputDynamics.sumCurrents(SimulationOptions.ContactNodes);
    OutputDynamics.networkResistance = Stimulus.Signal(ii,:)./OutputDynamics.sumCurrents(SimulationOptions.ContactNodes);
    %OutputDynamics.AverageLambda = avLambda;

    % ZK: also for local values:
    %OutputDynamics.lambda = lambda_vals;
    %OutputDynamics.storevoltage = voltage_vals;
    
end