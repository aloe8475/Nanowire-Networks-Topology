% 22/04/19
% Network Analysis Code v1.0
% --------------------------
dbstop if error
load_data_question=lower(input('Load network data, Analysis Data Only or None? N - None, D - Network Data, A - Analysis Data\n','s'));

if load_data_question=='d'
    clearvars -except load_data_question
    close all
    %load network data
    [network, network_load, simulations,sim_loaded]= load_data();
elseif load_data_question=='a'
    clear LDA_Analysis
    close all
    %Load previous LDA analysis data
    networkNum=input('Which Network # do you want to test? 1=500nw 2=100nw \n'); %% CHANGE WHICH SIMULATION YOU WANT TO TEST HERE.
    simNum=input('Which Simulation # do you want to test? \n'); %% CHANGE WHICH SIMULATION YOU WANT TO TEST HERE.
    [LDA_Analysis(simNum)] = load_LDA_data();
end

%%
%---
if load_data_question~='a'
    networkNum=input('Which Network # do you want to train? 1=500nw 2=100nw \n');
    if networkNum==1
        simNum=input('Which Simulation # do you want to train? (1-5) \n'); %% CHANGE WHICH SIMULATION YOU WANT TO TEST HERE.
    elseif networkNum==2
        simNum=6;
    end
end

fprintf(['Simulation: ' network(networkNum).Name simulations(simNum).Name ' selected \n\n']);

i = 1;
while i == 1
    analysis_type=lower(input('Which analysis would you like to perform? G - graph, L - LDA, N - none \n','s'));
    
    %% Graph Analysis
    if analysis_type=='g'
        %Call graph analysis function
        Graph=graph_analysis(network(networkNum),network_load,simulations,simNum);
        %Save Graph analysis data
        save_state=lower(input('Would you like to save the Graph Analysis? y or n \n','s'));
        if save_state=='y'
            save_graph(Graph,network(networkNum),network_load);
        end
        i=i+1;
    elseif analysis_type=='l'
        %% LDA Analysis
        %THIS IS WHERE YOU CHANGE THE VARIABLES FOR LDA
        drain1=[full(simulations(simNum).Data.IDrain1)]; %full(simulations(2).Data.IDrain1)];
        drain2=[full(simulations(simNum).Data.IDrain2)]; %full(simulations(2).Data.IDrain2)];
        %Current source
        source1=[full(simulations(simNum).Data.ISource1)]; %full(simulations(2).Data.ISource1)];
        source2=[full(simulations(simNum).Data.ISource2)]; %full(simulations(2).Data.ISource2)];
        %Voltage source
        source1Voltage=full(simulations(simNum).Data.VSource1);
        source2Voltage=full(simulations(simNum).Data.VSource2);
        Input = [source1 source2];
        Output = [drain1 drain2]; %OUTPUT is column (variable) x row (observations)
        Target = source1Voltage(:,1)>0.001;% uncomment for current: [source1(:,1)> 1.0e-04 *0.001 & source2<0];%TARGET is the classifier we expect
        
        %visualise drain1 and drain2
        figure
        plot(drain1); hold on
        plot(drain2);
        
        % Calculate linear discriminant coefficients - finding a line that
        % differentiate Output & Target
        LDA_Analysis(simNum).W = LDA(Output,Target);
        
        % % Calulcate linear scores for training data - what are the loadings for
        % that line - how would you get to 'x' from that line
        LDA_Analysis(simNum).L = [ones(size(Output,1),1) Output] * LDA_Analysis(simNum).W';
        %
        % % Calculate class probabilities
        LDA_Analysis(simNum).P = exp(LDA_Analysis(simNum).L) ./ repmat(sum(exp(LDA_Analysis(simNum).L),2),[1 2]);
        
        %Save all variables into struct
        LDA_Analysis(simNum).drain1=drain1;
        LDA_Analysis(simNum).drain2=drain2;
        LDA_Analysis(simNum).source1=source1;
        LDA_Analysis(simNum).source2=source2;
        LDA_Analysis(simNum).Output=Output;
        LDA_Analysis(simNum).Input=Input;
        LDA_Analysis(simNum).Target=Target;
        
        %% normalise LDA data:
        LDA_Analysis(simNum).normalisedOutput=LDA_normalise(Output);
        LDA_Analysis(simNum).normalisedInput=LDA_normalise(Input);
        
        LDA_Analysis(simNum).normalisedW=LDA(LDA_Analysis(simNum).normalisedOutput,Target);
        LDA_Analysis(simNum).normalisedL=[ones(size(LDA_Analysis(simNum).normalisedOutput,1),1) LDA_Analysis(simNum).normalisedOutput] * LDA_Analysis(simNum).normalisedW';
        LDA_Analysis(simNum).normalisedP=exp(LDA_Analysis(simNum).normalisedL) ./ repmat(sum(exp(LDA_Analysis(simNum).normalisedL),2),[1 2]);
        LDA_Analysis(simNum).TypeOfData='Training';
        clear drain1 drain2 source1 source2 Input Target Output
        i=i+1;%get out of while loop this loop finishes
        %% Saving LDA
        save_state=lower(input('Would you like to save the LDA Analysis? y or n \n','s'));
        if save_state=='y'
            save_LDA(LDA_Analysis(simNum),network(networkNum),network_load,sim_loaded);
            i=i+1; %get out of while loop this loop finishes
        end
        %% Plotting LDA
    end
    if analysis_type=='l' || analysis_type=='n'
        plot_state=lower(input('Would you like to plot LDA Analysis? y or n \n','s'));
        if plot_state=='y'
            plot_LDA(LDA_Analysis(simNum),simNum,network(networkNum).Name);
        end
        
        %% Apply LDA Training to testing data (different simulation)
        apply_LDA=lower(input('Would you like to apply the loaded LDA analysis to another network or simulation? y or n \n','s'));
        if apply_LDA=='y'
            networkNum=input('Which Network # do you want to test? 1=500nw 2=100nw \n');
            if networkNum==1
                simulationChoice=input('Which simulation # would you like to apply LDA to? (1-5) \n');
            elseif networkNum==2
                simulationChoice=6;
            end
            LDA_Analysis(simulationChoice).Output=[full(simulations(simulationChoice).Data.IDrain1) full(simulations(simulationChoice).Data.IDrain2)];
            LDA_Analysis(simulationChoice).Input=[full(simulations(simulationChoice).Data.ISource1) full(simulations(simulationChoice).Data.ISource2)];
            LDA_Analysis(simulationChoice).Target=[full(simulations(simulationChoice).Data.VSource1 >0.001)];
            LDA_Analysis(simulationChoice).TypeOfData='Testing';
            [LDA_Analysis(simulationChoice).appliedP, LDA_Analysis(simulationChoice).appliedL, LDA_Analysis(simulationChoice).normalisedOutput, LDA_Analysis(simulationChoice).normalisedInput]=LDA_Apply(LDA_Analysis(simNum).normalisedW,LDA_Analysis(simulationChoice).Output, LDA_Analysis(simulationChoice).Input);
            plot_state2=lower(input('Would you like to plot the applied LDA results? y or n \n','s'));
            if plot_state2=='y'
                plot_LDA_Applied(LDA_Analysis(simulationChoice),simNum,simulationChoice,network(networkNum).Name);
            end
        end
        i=i+1; %get out of while loop this loop finishes
    elseif analysis_type~='l' &&  analysis_type~='g' && analysis_type~='n'
        fprintf('Please type either G, L or N only');
    end
    
    if analysis_type=='g' || analysis_type=='n'
        plot_state=lower(input('Would you like to plot Graph Analysis? y or n \n','s'));
        if plot_state=='y'
            plot_graph(Graph,network(networkNum),simulations,sim_loaded,simNum)
        end
        i=i+1;
    end
    
    %% Insert Further Analysis Below
    % -------------------------------
    % -------------------------------
end

%--------------------------------------------------------------------------


%% FUNCTIONS
function plot_LDA(LDA_Analysis, simNum, networkName)
save_directory='D:\alon_\Research\POSTGRAD\PhD\CODE\Data\Figures\LDA\LDA Training\';
%plot LDA
LDAf=figure;
sgtitle('Electrode LDA Classification Training');
subplot(4,1,1)
plot(LDA_Analysis.Input)
title('Source Electrodes');
ylabel('Current (A)')
xlabel('Timestamp (0.01sec)')
subplot(4,1,2)
plot(LDA_Analysis.Output)
title('Drain Electrodes');
ylabel('Current (A)')
xlabel('Timestamp (0.01sec)')
subplot(4,1,3)
plot(LDA_Analysis.P(:,2))
title('LDA class probability (P)');
xlabel('Timestamp (0.01sec)')
ylabel('Probability');
subplot(4,1,4)
plot(LDA_Analysis.Target)
title('LDA Target Classification (class labels) - Source1 > 0V')
xlabel('Timestamp(0.01sec)');
ylabel('Voltage > 0');

%plot Log LDA
LDAff=figure;
sgtitle('Log-Scale LDA Electrode Classification Training');
subplot(4,1,1)
semilogy(LDA_Analysis.Input) %log y axis plot
title('Source Electrodes');
ylabel('Log Current (A)')
xlabel('Timestamp (0.01sec)')
subplot(4,1,2)
semilogy(LDA_Analysis.Output) %log y axis plot
title('Drain Electrodes');
ylabel('Log Current (A)')
xlabel('Timestamp (0.01sec)')
subplot(4,1,3)
plot(LDA_Analysis.P(:,2))
title('LDA class probability (P)');
xlabel('Timestamp (0.01sec)')
ylabel('Probability');
subplot(4,1,4)
plot(LDA_Analysis.Target)
title('LDA Target Classification (class labels) - Source1 > 0V')
xlabel('Timestamp(0.01sec)');
ylabel('Voltage > 0');

networkName(regexp(networkName,'[/:]'))=[]; %remove '/' character because it gives us saving problems

%Save
%Note; change the date of the simulation in the name if using a different
%set of simulations
saveas(LDAf,[save_directory num2str(networkName) '_6MaySim_LDA_Classification_Training_Simulation' num2str(simNum)],'jpg'); % CHANGE DATE OF SIMULATION
saveas(LDAf,[save_directory num2str(networkName) '_6MaySim_LDA_Classification_Training_Simulation' num2str(simNum)],'eps');
saveas(LDAff,[save_directory num2str(networkName) '_6MaySim_logLDA_Classification_Training_Simulation' num2str(simNum)],'jpg');
saveas(LDAff,[save_directory num2str(networkName) '_6MaySim_logLDA_Classification_Training_Simulation' num2str(simNum)],'eps');

end


function plot_LDA_Applied(LDA_Analysis,simNum,simulationChoice,networkName)
save_directory='D:\alon_\Research\POSTGRAD\PhD\CODE\Data\Figures\LDA\LDA Testing\';
appliedF=figure;
sgtitle('Electrode LDA Classification Testing');
subplot(4,1,1)
plot(LDA_Analysis.normalisedInput);
title('Normalised Source Electrodes');
ylabel('Current (A)')
xlabel('Timestamp (0.01sec)')
subplot(4,1,2)
plot(LDA_Analysis.normalisedOutput);
title('Normalised Drain Electrodes');
ylabel('Current (A)')
xlabel('Timestamp (0.01sec)')
subplot(4,1,3)
plot(LDA_Analysis.appliedP(:,2));
title('Classification Probability');
ylabel('Probability')
xlabel('Timestamp (0.01sec)')
subplot(4,1,4)
plot(LDA_Analysis.Target(:,1));
title('LDA Target Classification (class labels) - Source1 > 0V');
ylabel('Voltage > 0')
xlabel('Timestamp (0.01sec)')

networkName(regexp(networkName,'[/:]'))=[]; %remove '/' character because it gives us saving problems

saveas(appliedF,[save_directory num2str(networkName) '_6MaySim_LDA_Classification_TrainingSim_' num2str(simNum) '_TestingSim' num2str(simulationChoice)],'jpg');
saveas(appliedF,[save_directory num2str(networkName) '_6MaySim_LDA_Classification_TrainingSim_' num2str(simNum) '_TestingSim' num2str(simulationChoice)],'eps');
end

function [LDA_Analysis] = load_LDA_data()
cd('D:\alon_\Research\POSTGRAD\PhD\CODE\Data\LDA Analysis (Mac)');
waitfor(msgbox('Select the LDA Analysis saved data'));
[FileName,PathName] = uigetfile('*.mat','Select the LDA Analysis saved data');
f=fullfile(PathName,FileName);
load(f);
cd('D:\alon_\Research\POSTGRAD\PhD\CODE\Analysis');
end

function [network, network_load, simulations, sim_loaded] = load_data()

%% Load Data
%Ask to load Zdenka or Adrian:
network_load=lower(input('Which Network do you want to analyse? Z - Zdenka, A - Adrian \n','s'));

if strcmp(network_load,'a')
    %Get current network - Adrian
    [network,sim_loaded]=Load_Adrian_Code();
    %unpack simulation data into simulation variable
    if sim_loaded==1
        for i = 1:length(network.Simulations)
            simulations(i)=network.Simulations(i);
        end
    else
        simulations=network.Simulations;
    end
    cd('D:\alon_\Research\POSTGRAD\PhD\CODE\Analysis');
elseif strcmp(network_load,'z')
    %Get network - Zdenka:
    % D:\alon_\Research\POSTGRAD\PhD\CODE\Zdenka's Code\atomic-switch-network-1.3-beta\asn\connectivity\connectivity_data
    network=Load_Zdenka_Code();
    cd('D:\alon_\Research\POSTGRAD\PhD\CODE\Analysis');
end
end
function Graph=graph_analysis(network,network_load,simulations,simNum)
%% Mac's Analysis: (Graph)
if strcmp(network_load,'z')%Zdenka Code:
    net_mat=network.adj_matrix; %symmetrical matrix
elseif strcmp(network_load,'a') %adrian code
    IndexTime=input('What Timestamp do you want to analyse? 1-2000 \n'); %CHOOSE TIMESTAMP
    %this gives an adj matrix for the network used for a chosen simulation at a specific timestamp
    a= full(simulations(simNum).Data.Rmat{IndexTime});
    a(a==5000)=1;
    a(a==5000000)=0;  %binarising the resistance
    net_mat=a;
    %    net_mat=full(simulations(simNum).Data.AdjMat{IndexTime}); %this gives an adj matrix for the network used for a chosen simulation at a specific timestamp
end

%Global efficiency --> 1/characteristic path length, averaged over the whole network. An estimate of how integrated the network is.
Graph.GE = efficiency_bin(net_mat,0);

%Topological features will allow us to understand whether it was good or
%bad for classification.

%Local efficiency --> 1/characteristic path length, at each node
Graph.LE = efficiency_bin(net_mat,1);

%communicability --> an estimate of the ease with which each pair of nodes can connect in the network
Graph.COMM = expm(net_mat);

%modularity --> an estimate of how segregated the network is
[Graph.Ci,Graph.Q] = community_louvain(net_mat,1);
%The Ci term is the module assignment for each node
%The Q term is the 'quality' of the partition --> how modular the network is.
% -- this should tell us how well we can classify

%degree --> a count of how many edges are connected to each node
Graph.DEG = degrees_und(net_mat);

%participation coefficient --> an estimate of how integrative a node is
Graph.P = participation_coef(net_mat,Graph.Ci);
%Ci from 'community_louvain.m'

%module degree z-score --> an estimate of how segregated a node is
Graph.MZ = module_degree_zscore(net_mat,Graph.Ci);
%Ci from 'community_louvain.m'

% CIRCUIT RANK -- measure of recurrent loops (feedback loops)
% look in Zdenka's code for this

%save network matrix to graph struct
Graph.network=net_mat;
Graph.IndexTime=IndexTime;
end

function plot_graph(Graph, network, simulations,sim_loaded,simNum)
save_directory='D:\alon_\Research\POSTGRAD\PhD\CODE\Data\Figures\Graph Analysis\';

IndexTime=Graph.IndexTime;
%visualise graph network:
threshold=Graph.DEG>1; %greater than 1 degree threshold
Graph.networkThreshold=Graph.network(threshold,threshold); %applying degree threshold
g=graph(Graph.networkThreshold);

f1=figure;
p1=plot(g);

%find electrodes:
node_indices=find(threshold==1); %find nodes with threshold == 1
for i=1:4
    if ~isempty(find(node_indices==simulations(simNum).Electrodes.PosIndex(i)))
        new_electrodes(i).PosIndex=find(node_indices==simulations(simNum).Electrodes.PosIndex(i));
        new_electrodes(i).Name=simulations(simNum).Electrodes.Name(i);
    end
end

highlightElec={new_electrodes.PosIndex};
highlightElec=cell2num(highlightElec);
%highlight electrodes on graph
if sim_loaded==1
    highlight(p1,highlightElec,'NodeColor','green','MarkerSize',5); %change simulation number
    labelnode(p1,highlightElec,{new_electrodes(3).Name{1},new_electrodes(4).Name{1}}); %need to make this better - change 3:4 to a variable
end

f2=figure;
p2=plot(g);
%Closeness graph: (Unweighted)
ucc = centrality(g,'closeness');
p2.NodeCData=ucc;
colormap jet
colorbar
title(['Closeness Centrality Scores Timestamp ' num2str(IndexTime)])
if sim_loaded==1
    highlight(p2,highlightElec,'MarkerSize',7); %change simulation number
    labelnode(p2,highlightElec,{new_electrodes(3).Name{1},new_electrodes(4).Name{1}}); %need to make this better - change 3:4 to a variable

end

%Node size graph:
f3=figure;
p3=plot(g);
deg_ranks=Graph.DEG(threshold);
edges = linspace(min(deg_ranks),max(deg_ranks),7);
bins = discretize(deg_ranks,edges);
p3.MarkerSize = bins;
p3.NodeColor='r';
title(['Degree Size Timestamp ' num2str(IndexTime)]);
if sim_loaded==1
    highlight(p3,highlightElec,'NodeColor','green'); %change simulation number
    labelnode(p3,highlightElec,{new_electrodes(3).Name{1},new_electrodes(4).Name{1}}); %need to make this better - change 3:4 to a variable
end
text(-5,-6.2,'Min Degrees - 1 (small dot) | Max Degrees - 44 (large dot)');

%Both combined:
f4=figure;
p4=plot(g);
p4.MarkerSize = bins;
p4.NodeCData=ucc;
colormap jet
colorbar
labelnode(p3,highlightElec,{new_electrodes(3).Name{1},new_electrodes(4).Name{1}}); %need to make this better - change 3:4 to a variable
title(['Centrality Score and Degree Size Timestamp ' num2str(IndexTime)]);
text(-5.5,-6.2,'Min Degrees = 1 (small dot) | Max Degrees = 44 (large dot)');

%Histogram of degree distribution:
f5=figure;
h1=histogram(Graph.DEG(threshold));
title(['Distribution of Connectivity of Nodes Timestamp ' num2str(IndexTime)]);
xlabel('Number of Connections (Degrees)');
ylabel('Frequency');

%Cluster Analysis:
f6=figure;
p6=plot(g);
p6.MarkerSize = 4;
p6.NodeCData=Graph.Ci(threshold);
labelnode(p3,highlightElec,{new_electrodes(3).Name{1},new_electrodes(4).Name{1}}); %need to make this better - change 3:4 to a variable
colormap hsv(6) %change number of colors here if there are more/less than 6 clusters
title(['Cluster Analysis ' num2str(IndexTime)]);

%Participant Coefficient Analysis:
f7=figure;
p7=plot(g);
p7.NodeCData=Graph.Ci(threshold);
p_ranks=Graph.P(threshold);
edges2 = linspace(min(p_ranks),max(p_ranks),7);
bins2 = discretize(p_ranks,edges2);
p7.MarkerSize=bins2;
colormap hsv(6)
labelnode(p3,highlightElec,{new_electrodes(3).Name{1},new_electrodes(4).Name{1}}); %need to make this better - change 3:4 to a variable
text(-5.5,-6.2,['Min P Coeff = 0 (small dot) | Max P Coeff = ' num2str(max(Graph.P(threshold))) ' (large dot)']);
title(['Participant Coefficient Analysis Timestamp ' num2str(IndexTime)]);

%Module Degree Z-Score:
f8=figure;
p8=plot(g);
p8.NodeCData=Graph.Ci(threshold);
mod_ranks=Graph.MZ(threshold);
edges3 = linspace(min(mod_ranks),max(mod_ranks),7);
bins3 = discretize(mod_ranks,edges3);
p8.MarkerSize=bins3;
colormap hsv(6)
labelnode(p3,highlightElec,{new_electrodes(3).Name{1},new_electrodes(4).Name{1}}); %need to make this better - change 3:4 to a variable
text(-6,-6.2,['Min MZ Coeff = ' num2str(min(Graph.MZ(threshold))) ' (small dot) | Max P Coeff = ' num2str(max(Graph.MZ)) ' (large dot)']);
title(['Module Degree z-Score Analysis Timestamp ' num2str(IndexTime)]);


% Communicability at different times:
Adj=(simulations(simNum).Data.AdjMat{IndexTime});%convert 498x498 matrix to EdgeCData ~(1x6065)
Adj=Adj(threshold,threshold);
[j,i,~]=find(tril(Adj));
wd=zeros(1,length(j));
Graph.COMM=Graph.COMM(threshold,threshold);
for k=1:length(j)
    wd(k)=Graph.COMM(i(k),j(k));
end

Adj2=(simulations(simNum).Data.AdjMat{IndexTime});%convert 498x498 matrix to EdgeCData ~(1x6065)
Adj2=Adj2(threshold,threshold);
[j,i,~]=find(tril(Adj2));
wd2=zeros(1,length(j));

for k=1:length(j)
    wd2(k)=Graph.networkThreshold(i(k),j(k));
end

wd3=wd(logical(wd2));
f9=figure;
p9=plot(g);
p9.EdgeCData=wd3;%log10(wd);
p9.MarkerSize=2;
colormap jet
colorbar
labelnode(p3,highlightElec,{new_electrodes(3).Name{1},new_electrodes(4).Name{1}}); %need to make this better - change 3:4 to a variable
title(['Communicability Analysis Timestamp ' num2str(IndexTime) ' (log10)']);

%Save
network.Name(regexp(network.Name,'[/:]'))=[]; %remove '/' character because it gives us saving problems

saveas(f1,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_Graph_Network_Timestamp' num2str(IndexTime)],'jpg');
saveas(f1,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_Graph_Network_Timestamp' num2str(IndexTime)],'eps');
saveas(f2,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_Graph_Closeness_Centrality_Timestamp' num2str(IndexTime)],'jpg');
saveas(f2,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_Graph_Closeness_Centrality_Timestamp' num2str(IndexTime)],'eps');
saveas(f3,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_Graph_Degree_Size_Timestamp' num2str(IndexTime)],'jpg');
saveas(f3,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_Graph_Degree_Size_Timestamp' num2str(IndexTime)],'eps');
saveas(f4,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_Graph_Degree_&_Closeness_Timestamp' num2str(IndexTime)],'jpg');
saveas(f4,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_Graph_Degree_&_Closeness_Timestamp' num2str(IndexTime)],'eps');
saveas(f5,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_Graph_Degree_Distribution_Timestamp' num2str(IndexTime)],'jpg');
saveas(f5,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_Graph_Degree_Distribution_Timestamp' num2str(IndexTime)],'eps');
saveas(f6,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_Graph_Cluster_Analysis_Timestamp' num2str(IndexTime)],'jpg');
saveas(f6,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_Graph_Cluster_Analysis_Timestamp' num2str(IndexTime)],'eps');
saveas(f7,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_Graph_Participant_Coefficient_Analysis_Timestamp' num2str(IndexTime)],'jpg');
saveas(f7,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_Graph_Participant_Coefficient_Analysis_Timestamp' num2str(IndexTime)],'eps');
saveas(f8,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_Graph_Module_Degree_Z_Score_Analysis_Timestamp' num2str(IndexTime)],'jpg');
saveas(f8,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_Graph_Module_Degree_Z_Score_Analysis_Timestamp' num2str(IndexTime)],'eps');
saveas(f9,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_Graph_Communicability_Analysis_Timestamp' num2str(IndexTime)],'jpg');
saveas(f9,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_Graph_Communicability_Analysis_Timestamp' num2str(IndexTime)],'eps');
end

function save_graph(Graph,network,network_load)
save_directory='D:\alon_\Research\POSTGRAD\PhD\CODE\Data\Graph Analysis (Mac)\';
if strcmp(network_load,'z')%Zdenka Code:
    save([save_directory 'Zdenka_' num2str(network.number_of_wires) 'nw_Graph_Analysis_' date],'Graph');
elseif strcmp(network_load,'a') %adrian code
    network.Name(regexp(network.Name,'[/:]'))=[]; %remove '/' character because it gives us saving problems
    save([save_directory 'Adrian_' num2str(network.Name) 'Graph_Analysis_' date],'Graph');
end
end

function save_LDA(LDA_Analysis,network,network_load,sim_loaded)
save_directory='D:\alon_\Research\POSTGRAD\PhD\CODE\Data\LDA Analysis (Mac)\';
if strcmp(network_load,'z')%Zdenka Code:
    save([save_directory 'Zdenka_' num2str(network.number_of_wires) 'nw_LDA_Analysis_' date],'LDA_Analysis');
elseif strcmp(network_load,'a') %adrian code
    network.Name(regexp(network.Name,'[/:]'))=[]; %remove '/' character because it gives us saving problems
    save([save_directory 'Adrian_' num2str(network.Name) 'LDA_Analysis_' date],'LDA_Analysis','sim_loaded');
end
end