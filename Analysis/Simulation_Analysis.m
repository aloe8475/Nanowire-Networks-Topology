

function Simulation_Analysis() 
%load Adrian's code
network_Adrian=Load_Adrian_Code();

%unpack simulation data into simulation variable
for i = 1:length(network_Adrian.Simulations)
simulations(i)=network_Adrian.Simulations{i};
end 
%% Simulation Data Analysis:
%Simulation data consists of Data Tables with the following columns:

%AdjMat, Gmat, Wmat, Rmat, SumRule, Currents, Voltages, ISource1, VSource1,
%IDrain1, VDrain1, ISource2, VSource2, IDrain2, VDrain2, VoltDif, time. 

%Note: for more electrodes, there will be more columns with the electrode
%outputs.

%Check data:
Data.check=head(simulations(1).Data); %show first eight rows of the Data Table (peak at the data)

%Statistical Summary:
Data.summary=summary(simulations(1).Data);

end 