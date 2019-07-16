    %% Across Network Exploration:
    %-------------------------------------------------------------------------

    % This script loads Graph Theory analyses performed in Network_Explore for
    % all sizes of networks we have simulated. It then creates random and
    % ordered networks of the same size (bootstrapped x100), and compares Graph Theoretical
    % concepts between our networks and the random/ordered networks.

    % Author: Alon Loeffler
    % Date: 19/06/2019
    % Version: 1.1

    % Changelog:
    % 26/06/2019 - Added Participant coefficient, Small World Propensity and
    % Betweenness centrality measures.
    %-------------------------------------------------------------------------

    clear all
    close all

    computer=getenv('computername');
    switch computer
        case 'W4PT80T2' %if on desktop at uni - Alon
            explore_location='C:\Users\aloe8475\Documents\PhD\GitHub\CODE\Data\Explore Analysis\';
        case '' %if on linux
            explore_location='/suphys/aloe8475/Documents/CODE/Data/Explore Analysis/';
        case 'LAPTOP-S1BV3HR7'
            explore_location='D:\alon_\Research\POSTGRAD\PhD\CODE\Data\Explore Analysis\';
            %case '' %--- Add other computer paths (e.g. Mike)
    end

    %Load three explore analyses:
    e100=load([explore_location 'Adrian_Net_Sx20_NoW100_0325-2019_112338__Sim_1_SourceElectrode_6_DrainElectrode_76_Exploration_Analysis_ Timestamp_400_26-Jun-2019.mat']);
    e500=load([explore_location 'Adrian_Net_Sx20_NoW500_0330-2019_111659__Sim_1_SourceElectrode_18_DrainElectrode_430_Exploration_Analysis_ Timestamp_400_26-Jun-2019.mat']);
    e1000=load([explore_location 'Adrian_Net_Sx20_NoW1000_0606-2019_113353__Sim_1_SourceElectrode_32_DrainElectrode_1000_Exploration_Analysis_ Timestamp_400_26-Jun-2019.mat']);
    e2000=load([explore_location 'Adrian_Net_Sx20_NoW2000_0618-2019_125103__Sim_1_SourceElectrode_158_DrainElectrode_1820_Exploration_Analysis_ Timestamp_400_26-Jun-2019.mat']);

    %% Example AI Graph Analysis (Recurrent Neural Network)
    
    %Create Sample RNN: - NEED TO TALK TO MAC ABOUT THIS
    %AI.AdjMat=zeros([500 500]);
    
    %% C-Elegans:
    cElegans.AvgMZ=0.62;% - Guimera & Nunes, 2005
    
    %% Human Graph Analysis
   %Cluster Coeff & Path Length
    human.GlobalClust=0.53; %Taken from (Achard et al., 2006) 
    human.AvgPath=2.49; %Taken from (Achard et al., 2006)
    %Participation Coefficient & Module z-Score
    human.PLocalHubs=0.2; %High PCoeff - approximates from (Power et al., 2013)
    human.PConnectorHubs=0.9; %Low PCoeff
    human.AvgP=0.502; %Bertolero, Yeo & D'Esposito - 2015
    human.MZHubs=median([2.5:0.5:10]);%[2.5 10] % within module degree probability less than 0.01 is analogous to having a z-score above 2.5 (Guimera and Amaral, 2005, Joyce et al., 2010)
    human.MZNonHubs=median([-2:0.5:2.5]);% (Guimera and Amaral, 2005, Joyce et al., 2010)
    human.AvgMZ=0.0001; %Bertolero, Yeo & D'Esposito - 2015
    %Small world Prop
    
    %Complexity
    
    %Communicability
    
    %Betweenness Centrality
    
    %Circuit Rank

    %% Random and Ordered Graph Analysis

    createNewRand=lower(input('Would you like to create new Random and Ordered graphs? (Note this will take 4+ Hours) \n','s'));
    loadPath=explore_location;
    if createNewRand=='n'
        if exist([loadPath 'Ordered_Graphs_500nw.mat'], 'file') == 2 && exist([loadPath 'Random_Graphs_500nw.mat'],'file') == 2 %2 because .mat file
            load([loadPath 'Ordered_Graphs_500nw.mat']);
            load([loadPath 'Random_Graphs_500nw.mat']);
            fprintf('Loading Ordered and Random Graphs (500nodes)');
        else
            fprintf('Ordered and Random Graphs have not been created yet \n');
            fprintf('Creating New Graphs \n');
            [random, ordered, random100, ordered100]=createRandom_Ordered_Graphs(e500);
        end
    else
        if exist([loadPath 'Ordered_Graphs_500nw.mat'], 'file') == 2 && exist([loadPath 'Random_Graphs_500nw.mat'],'file') == 2 %2 because .mat file
            overwrite=lower(input('Graphs (500nodes) have already been created, would you like to overwrite? \n','s'));
            if overwrite =='n'
                fprintf('Loading Ordered and Random Graphs (500nodes) \n');
                load([loadPath 'Ordered_Graphs_500nw.mat']);
                load([loadPath 'Random_Graphs_500nw.mat']);
            else
                fprintf('Creating New Graphs \n');
                [random, ordered, random100, ordered100]=createRandom_Ordered_Graphs(e500);
                fprintf('New Graphs Created \n');
            end
        else
            fprintf('Creating New Graphs \n');
            [random, ordered, random100, ordered100]=createRandom_Ordered_Graphs(e500);
            fprintf('New Graphs Created \n');
        end
    end

    fprintf('Graphs Loaded... \n Creating Plots... \n');


    %% Random & Ordered
    %Number of Nodes
    for i =1:100
        random(i).numNodes=height(random(i).Graph.Nodes);
        ordered(i).numNodes=height(ordered(i).Graph.Nodes);
        random(i).numEdges=height(random(i).Graph.Edges);
        ordered(i).numEdges=height(ordered(i).Graph.Edges);
    end 
    random100.AvgNodes=mean([random(:).numNodes]);
    ordered100.AvgNodes=mean([ordered(:).numNodes]);
    
    %Number of Edges
    random100.AvgEdges=mean([random(:).numEdges]);
    ordered100.AvgEdges=mean([ordered(:).numEdges]);
    
    % participation coefficient (mean)
    % Avg & std PCoeff per Node (across 100 bootstraps) and across nodes:
    random100.AvgAvgPCoeff=mean(mean([random(:).P]),2);
    random100.StdAvgPCoeff=std(mean([random(:).P]),[],2);
    ordered100.AvgAvgPCoeff=mean(mean([ordered(:).P]),2);
    ordered100.StdAvgPCoeff=std(mean([ordered(:).P]),[],2);
    
    %Small World Prop
    random100.AvgSmallWorldProp=mean([random(:).SmallWorldProp],2);
    random100.StdSmallWorldProp=std([random(:).SmallWorldProp],[],2);
    ordered100.AvgSmallWorldProp=mean([ordered(:).SmallWorldProp],2);
    ordered100.StdSmallWorldProp=std([ordered(:).SmallWorldProp],[],2);
    
    % communicability
    
    % complexity (from the Sporns, Tononi and Edelman paper I sent through the other day -- still waiting on code from Olaf, but will send through when I get it).
    
    % betweeness_centrality

    %Average Normalised Betweenness (BC/[(N-1)(N-2)])
%     random100.AvgNormBC=mean([random(:).normBC]);
%     random100.StdNormBC=std([random(:).normBC]);
%     ordered100.AvgNormBC=mean([ordered(:).normBC]);
%     ordered100.StdNormBC=std([ordered(:).normBC]);
   
    %% AgNW
    %Circuit Rank
    AgNW.CircuitRank=[e100.Explore.GraphTheory.CircuitRank e500.Explore.GraphTheory.CircuitRank e1000.Explore.GraphTheory.CircuitRank e2000.Explore.GraphTheory.CircuitRank];
    AgNW.GlobalClust=[e100.Explore.GraphTheory.GlobalClust, e500.Explore.GraphTheory.GlobalClust, e1000.Explore.GraphTheory.GlobalClust e2000.Explore.GraphTheory.GlobalClust];
    AgNW.AvgPath=[e100.Explore.GraphTheory.AvgPath, e500.Explore.GraphTheory.AvgPath, e1000.Explore.GraphTheory.AvgPath e2000.Explore.GraphTheory.AvgPath];
    AgNW.SmallWorldProp=[e100.Explore.GraphTheory.SmallWorldProp, e500.Explore.GraphTheory.SmallWorldProp, e1000.Explore.GraphTheory.SmallWorldProp e2000.Explore.GraphTheory.SmallWorldProp];

    %% Plot:
    % Small World Analysis
    x=[random100.AvgGlobalClust human.GlobalClust ordered100.AvgGlobalClust e100.Explore.GraphTheory.GlobalClust e500.Explore.GraphTheory.GlobalClust e1000.Explore.GraphTheory.GlobalClust e2000.Explore.GraphTheory.GlobalClust];
    y=[random100.AvgPath human.AvgPath ordered100.AvgPath e100.Explore.GraphTheory.AvgPath, e500.Explore.GraphTheory.AvgPath, e1000.Explore.GraphTheory.AvgPath e2000.Explore.GraphTheory.AvgPath];
    f=figure;
    p=gscatter(x,y);
    hold on
    e=errorbar(x(1), y(1),random100.StdPath);
    e2=errorbar(x(1), y(1),random100.StdGlobalClust,'horizontal');
    errorbar(x(3), y(3),ordered100.StdPath);
    errorbar(x(3), y(3),ordered100.StdGlobalClust,'horizontal');
    % xlim([0.05 0.6])
    % ylim([2 16])
    text(x,y,{'500node Random Nw','Human Nw','500node Ordered Nw','100nw','500nw','1000nw','2000nw'},'VerticalAlignment','bottom','HorizontalAlignment','left')
    xlabel('Global Clustering Coefficient');
    ylabel('Global Mean Path Length');
    p.MarkerEdgeColor='b';
    p(:,1).MarkerEdgeColor='r';
    p.LineWidth=1.5;

    %Small World Log
    f1=figure;
    logx=log10(x);
    logy=log10(y);
    p1=gscatter(x,logy);
    hold on
    hold on
    e=errorbar(x(1), logy(1),random100.StdPath);
    e2=errorbar(x(1), logy(1),random100.StdGlobalClust);
    errorbar(x(3), logy(3),ordered100.StdPath);
    errorbar(x(3), logy(3),ordered100.StdGlobalClust);
    % xlim([0.05 0.6])
    % ylim([2 16])
    ylim([0,max(logy)]);
    text(x,logy,{'500node Random Nw','Human Nw','500node Ordered Nw','100nw','500nw','1000nw','2000nw'},'VerticalAlignment','bottom','HorizontalAlignment','left')
    xlabel('Global Clustering Coefficient');
    ylabel('Log10 Global Mean Path Length');
    p1.MarkerEdgeColor='b';
    p1(:,1).MarkerEdgeColor='r';
    p1.LineWidth=1.5;

    %Small World Prop
    x=[random100.AvgSmallWorldProp ordered100.AvgSmallWorldProp AgNW.SmallWorldProp];
    ff=figure;
    pp=bar(x);
    hold on
    e=errorbar(x(1), random100.StdSmallWorldProp);
    e2=errorbar(x(2),ordered100.StdSmallWorldProp);
    % xlim([0.05 0.6])
    % ylim([2 16])
    xticklabels({'500node Random Nw','500node Ordered Nw','100nw','500nw','1000nw','2000nw'});
    ylabel('Small World Prop');
    
    %Circuit Rank:
    circuitRank=[random(1).CircuitRank ordered(1).CircuitRank AgNW.CircuitRank];
    f2=figure;
    p2=bar(circuitRank);
    xticklabels({'500node Random Nw', '500node Ordered Nw', '100nw', '500nw', '1000nw','2000nw'});
    ylabel('Circuit Rank');



    hAx=gca;            % get a variable for the current axes handle
    hT=[];              % placeholder for text object handles
    for i=1:length(p2)  % iterate over number of bar objects
        hT=[hT text(p2(i).XData+p2(i).XOffset,p2(i).YData,num2str(p2(i).YData.'), ...
            'VerticalAlignment','bottom','horizontalalign','center')];
    end

    %Average Participation Coefficient & Module z-Score
    f3=figure;
    %High PCoeff = Hubs / Central areas (Power et al., 2013)
    PCoeff=[random100.AvgAvgPCoeff human.PLocalHubs human.PLocalHubs human.PConnectorHubs human.PConnectorHubs ordered100.AvgAvgPCoeff mean(e100.Explore.GraphTheory.P), mean(e500.Explore.GraphTheory.P), mean(e1000.Explore.GraphTheory.P) mean(e2000.Explore.GraphTheory.P)];
    MZ=[random100.AvgAvgMZ human.MZHubs human.MZNonHubs human.MZHubs human.MZNonHubs ordered100.AvgAvgMZ mean(e100.Explore.GraphTheory.MZ), mean(e500.Explore.GraphTheory.MZ), mean(e1000.Explore.GraphTheory.MZ) mean(e2000.Explore.GraphTheory.MZ)];
    p3=gscatter(PCoeff,MZ);
    text(PCoeff,MZ,{'500node Random Nw', 'Human Connector Local Provincial Hub','Human Local Peripheral Node','Human Connector Hub','Human Satellite Connector', '500node Ordered Nw', '100nw', '500nw', '1000nw','2000nw'});
    xlabel('Average Participant Coefficient Coefficient');
    ylabel('Average Module z-Score');
    p3.MarkerEdgeColor='b';
    p3(:,1).MarkerEdgeColor='r';
    p3.LineWidth=1.5;
    
    %500 Node Participation Coefficient & Module z-Score
    f4=figure;
    %High PCoeff = Hubs / Central areas (Power et al., 2013)
    PCoeff=[random100.AvgPCoeff ordered100.AvgPCoeff, e500.Explore.GraphTheory.P];
    MZ=[random100.AvgMZ  ordered100.AvgMZ, e500.Explore.GraphTheory.MZ];
    p3=gscatter(PCoeff,MZ);
    text(PCoeff,MZ,{'500node Random Nw', '500node Ordered Nw', '500nw'});
    xlabel('Participant Coefficient Coefficient');
    ylabel('Module z-Score');
    p3.MarkerEdgeColor='b';
    p3(:,1).MarkerEdgeColor='r';
    p3.LineWidth=1.5;
    
    %Communicability
    
    %Complexity
    
    %Betweenness Centrality
    

    %% OLD CODE
    %% Load three networks
    % n100=load('Net_Sx_20_NoW100_03_25-2019_11_23_38_.mat');
    % n500=load('Net_Sx_20_NoW500_03_30-2019_11_16_59_.mat');
    % n1000=load('Net_Sx_20_NoW1000_06_06-2019_11_33_53_.mat');
    %
    % %Load three simulations:
    % s100=load('Net_Sx_20_NoW100_03_25-2019_11_23_38_Zdenka_Constant_6SimsOnly_4_Sec_Vmax_1_20-May-2019.mat');
    % s500=load('Net_Sx_20_NoW500_03_30-2019_11_16_59_Zdenka_Constant_1SimsOnly_4_Sec2Electrodes_Vmax_0.25_06-Jun-2019.mat');
    % s1000=load('Net_Sx_20_NoW1000_06_06-2019_11_33_53_Zdenka_Constant_1SimsOnly_4_Sec2Electrodes_Vmax_0.25_06-Jun-2019.mat');
    %
    %
    % %Set up 100nw network
    % if exist('s100.SelSims','var') == 1
    %       temp= s100.SelSims;%save simulations from test network if it hasn't been changed
    % else
    %       temp = n100.network.Simulations;%save simulations from test network if it has been changed
    % end
    %
    % n100.network.Simulations=temp;
    %
    % %Set up 500nw Network
    % if exist('s500.SelSims','var') == 1
    %       temp= s500.SelSims;%save simulations from test network if it hasn't been changed
    % else
    %       temp = n500.network.Simulations;%save simulations from test network if it has been changed
    % end
    %
    % n500.network.Simulations=temp;
    %
    % %Set up 1000nw Network
    % if exist('s1000.SelSims','var') == 1
    %       temp= s1000.SelSims;%save simulations from test network if it hasn't been changed
    % else
    %       temp = n1000.network.Simulations;%save simulations from test network if it has been changed
    % end
    %
    % n1000.network.Simulations=temp;