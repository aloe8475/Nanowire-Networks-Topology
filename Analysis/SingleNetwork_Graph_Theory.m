% 06/10/19
% Network Analysis Code - This code performs graph and subgraph analysis for
% individual Explore files. 
% Use this instead of MultiNetwork_Graph_Theory.m when we need to explore
% different simulations of 1 network. Use together with
% MultiSims_Graph_Theory.m 
% --------------------------
dbstop if error
load_data_question=lower(input('Load Explore Data or None? N - None, D - Explore Data \n','s'));

if load_data_question=='d'
    clearvars -except load_data_question
    close all
    %load network data
    [Explore,ws]= load_data();
end
    
sizeNetwork=height(Explore.View.Nodes);

% DEGREE
fDeg=figure;
    h=histogram(Explore.GraphTheory.DEG);
    title([num2str(sizeNetwork) 'nw']);
    xlabel('Degree');
    ylabel('Frequency');
title('Degree Distribution Comparison of Nanowire Networks')


%% WATTS STROGATZ
% Small World Analysis

% loadPath='C:\Users\aloe8475\Documents\PhD\GitHub\CODE\Analysis\Watson Strogatz\';
% b=load([loadPath 'beta.mat']);
% watStr.b=b.beta;
% load([loadPath 'cc.mat']);
% watStr.cc=cc;
% load([loadPath 'pl.mat']);
% watStr.pl=pl;

fWatts=figure('Position',[0 0 1920 1080]);

%Watss-Strogatz Colors:
lightblue=rgb('cyan');
blue=rgb('deep blue');
lightred=rgb('rose pink');
red=rgb('burnt red');


rdif=[blue(1)-lightblue(1)];
r2=blue(1):(-rdif)/100:lightblue(1);
gdif=[blue(2)-lightblue(2)];
g2=blue(2):(-gdif)/100:lightblue(2);
bdif=[blue(3)-lightblue(3)];
b2=blue(3):(-bdif)/100:lightblue(3);

rdif=[red(1)-lightred(1)];
r3=red(1):(-rdif)/100:lightred(1);
gdif=[red(2)-lightred(2)];
g3=red(2):(-gdif)/100:lightred(2);
bdif=[red(3)-lightred(3)];
b3=red(3):(-bdif)/100:lightred(3);

clrs2={[r2; g2; b2]'};
clrs3={[r3; g3; b3]'};

% s(1)=scatter(watStr.pl,mean(watStr.cc),[],[r3;g3;b3;]');
hold on

% f=figure;

   
%     if plotN3 %if we only want 100 & 500 networks
        y=[randomOrdered100.AvgGlobalClust Explore.GraphTheory.GlobalClust];
        x=[randomOrdered100.AvgPath Explore.GraphTheory.AvgPath];
        for i=1:length(x)
            %           s(i+1)=scatter(nanmean(categories.Clust{i}),categories.PathLength{i,j},[],clrs2{1}(i,:));
            s(i)=scatter(x(:,i),y(:,i),[],clrs2{1}(i,:));
            hold on
            e=errorbar(x(1), y(1),Net(plotNet).randomOrdered100.StdPath);
            e2=errorbar(x(1), y(1),Net(plotNet).randomOrdered100.StdGlobalClust,'horizontal');
%             e3=errorbar(x(5), y(5),Net(plotNet).ordered100.StdPath);
%             e4=errorbar(x(5), y(5),Net(plotNet).ordered100.StdGlobalClust,'horizontal');
            e.Color=clrs2{1}(i,:);
            e2.Color=clrs2{1}(i,:);
%             e3.Color=clrs2{1}(i,:);
%             e4.Color=clrs2{1}(i,:);
        end
            text(x,y,{[num2str(sizeNetwork) 'node Watts-Strogats'] [num2str(sizeNetwork) 'nw']},'VerticalAlignment','bottom','HorizontalAlignment','left')

    ylabel('Global Clustering Coefficient');
    xlabel('Global Mean Path Length');
    p.MarkerEdgeColor='b';
    p(:,1).MarkerEdgeColor='r';
    p.LineWidth=1.5;


fprintf('Figure 1 Complete \n');

smallWorldPropOrderedRandom=[];


    %Circuit Rank:
    smallWorldPropOrderedRandom=[smallWorldPropOrderedRandom randomOrdered100.AvgSmallWorldProp];
%     smallWorldPropOrdered=[smallWorldPropOrdered Net(i).ordered100.AvgSmallWorldProp];
    randomLabel=[num2str(sizeNetwork) ' Watts-Strogatz Nw'];
%     orderedLabel{i}=[num2str(Net(i).sizeNetwork) ' Ordered Nw'];


%Small World Prop
f2=figure;
%smallWorldPropOrdered
x=[smallWorldPropOrderedRandom  Explore.GraphTheory.SmallWorldProp];
p2=bar(x);
hold on
    e=errorbar(x(1), randomOrdered100.StdSmallWorldProp);
%     e2=errorbar(x(2),Net(i).ordered100.StdSmallWorldProp);
    % xlim([0.05 0.6])
    % ylim([2 16])
    hold on
%orderedLabel
xticklabels([randomLabel,[num2str(sizeNetwork) 'nw']]);

set(gca, 'XTickLabelRotation', 45)

ylabel('Small World Prop');


fprintf('Figure 2 Complete \n');


f3=figure;
circuitRankRandom=[];
% circuitRankOrdered=[];

    %Circuit Rank:
    circuitRankRandom=[circuitRankRandom randomOrdered(1).CircuitRank];
%     circuitRankOrdered=[circuitRankOrdered Net(i).ordered(1).CircuitRank];


circuitRank=[circuitRankRandom];
p3=bar(circuitRank);
% orderedLabel
xticklabels([randomLabel, [num2str(sizeNetwork) 'nw']]);
set(gca, 'XTickLabelRotation', 45)

ylabel('Circuit Rank');
hAx=gca;            % get a variable for the current axes handle
hT=[];              % placeholder for text object handles
for i=1:length(p3)  % iterate over number of bar objects
    hT=[hT text(p3(i).XData+p3(i).XOffset,p3(i).YData,num2str(p3(i).YData.'), ...
        'VerticalAlignment','bottom','horizontalalign','center')];
end
hold on
%Average Participation Coefficient & Module z-Score
fprintf('Figure 3 Complete \n');

f4=figure;

PRandom=[];
% POrdered=[];
MZRandom=[];
% MZOrdered=[];
stdPRandom=[];
% stdPOrdered=[];
stdMZRandom=[];
% stdMZOrdered=[];

    %Circuit Rank:
    PRandom=[PRandom randomOrdered100.AvgAvgPCoeff];
%     POrdered=[POrdered Net(i).ordered100.AvgAvgPCoeff];
    stdPRandom=[stdPRandom randomOrdered100.StdAvgPCoeff];
%     stdPOrdered=[stdPOrdered Net(i).ordered100.StdAvgPCoeff];
    MZRandom=[MZRandom randomOrdered100.AvgAvgMZ];
%     MZOrdered=[MZOrdered Net(i).ordered100.AvgAvgMZ];
    stdMZRandom=[stdMZRandom randomOrdered100.StdAvgMZ];
%     stdMZOrdered=[stdMZOrdered Net(i).ordered100.StdAvgMZ];
% POrdered stdMZOrdered MZOrdered

PCoeff=[PRandom ,  mean(Explore.GraphTheory.P)];
stdPCoeff=[stdPRandom,  std(Explore.GraphTheory.P)];
stdMZ=[stdMZRandom std(Explore.GraphTheory.MZ)];
    MZ=[MZRandom mean(Explore.GraphTheory.MZ)];
    
p4=gscatter(PCoeff,MZ);
% e=errorbar(MZ, stdMZ);
% e2=errorbar(PCoeff,stdPCoeff);

%High PCoeff = Hubs / Central areas (Power et al., 2013)
%stdMZRandom
text(PCoeff,MZ,[randomLabel, [num2str(sizeNetwork) 'nw']],'NorthWest');
xlabel('Average Participant Coefficient Coefficient');
ylabel('Average Module z-Score');
p4.MarkerEdgeColor='b';
p4(:,1).MarkerEdgeColor='r';
p4.LineWidth=1.5;
hold on

fprintf('Figure 4 Complete \n');

%% Plot Guimera & Amaral rectangles: ALON TO DO 
f5=figure;
    guimera(network,Net,i); %change network here
fprintf('Figure 5 Complete \n');

% %Communicability
% f6=figure;
% COMM=[Net(plotNet).randomOrdered100.AvgCOMM(1) Net(plotNet).ordered100.AvgCOMM(1) AgNW.AvgCOMM];
% stdCOMM=[0 0 AgNW.StdCOMM];
% p6=bar(log10(COMM));
% hold on
% e=errorbar(log10(COMM), log10(stdCOMM)); %use log10 otherwise communicability is much too large to visualise
% e.LineStyle='none';
% % xlim([0.05 0.6])
% % ylim([2 16])
% xticklabels({[num2str(Net(plotNet).sizeNetwork) 'node Random Nw'],[num2str(Net(plotNet).sizeNetwork) 'node Ordered Nw'],'100nw','500nw','1000nw','2000nw'});
% ylabel('Log10 Communicability');
fprintf('Figure 6 Complete \n');

%Betweenness Centrality

f7=figure;
BCRandom=[];
BCOrdered=[];
BCRandomstd=[];
BCOrderedstd=[];
    %Circuit Rank:
    BCRandom=[BCRandom randomOrdered100.AvgBC];
%     BCOrdered=[BCOrdered Net(i).ordered100.AvgBC];
    BCRandomstd=[BCRandomstd randomOrdered100.StdBC];
%     BCOrderedstd=[BCOrderedstd Net(i).ordered100.StdBC];
    randomLabel{i}=[num2str(sizeNetwork) ' Nw'];
%     orderedLabel{i}=[num2str(Net(i).sizeNetwork) ' Ordered Nw'];


BC=[BCRandom BCOrdered ANN.avgBC cElegans.avgBC AgNW.AvgBC];
stdBC=[BCRandomstd BCOrderedstd ANN.stdBC cElegans.stdBC AgNW.StdBC];
p7=bar(BC);
hold on
e=errorbar(BC, stdBC);
e.LineStyle='none';
% xlim([0.05 0.6])
% ylim([2 16])

xticklabels([randomLabel{:},orderedLabel{:},'500node Artificial Neural Nw','C. Elegans Nw','100nw','500nw','1000nw','2000nw']);
set(gca, 'XTickLabelRotation', 45)
ylabel('Betweenness Centrality');
hold on


fprintf('Figure 7 Complete \n');



%% SAVE GRAPHS
set(f,'PaperPositionMode','auto');
set(f,'PaperOrientation','landscape');
set(f,'Position',[0 0 1920 1080]);
print(f,'-painters','-dpdf','-bestfit','-r600',[fig_dir 'Clustering Coefficient vs Path Length all networks.pdf']);

set(fDeg,'PaperPositionMode','auto');
set(fDeg,'PaperOrientation','landscape');
set(fDeg,'Position',[0 0 1920 1080]);
print(fDeg,'-painters','-dpdf','-bestfit','-r600',[fig_dir 'Degree Distribution all networks.pdf']);


set(f2,'PaperPositionMode','auto');
set(f2,'PaperOrientation','landscape');
set(f2,'Position',[0 0 1920 1080]);
print(f2,'-painters','-dpdf','-bestfit','-r600',[fig_dir 'Small World Prop all networks.pdf']);

set(f3,'PaperPositionMode','auto');
set(f3,'PaperOrientation','landscape');
set(f3,'Position',[0 0 1920 1080]);
print(f3,'-painters','-dpdf','-bestfit','-r600',[fig_dir 'Circuit Rank all networks.pdf']);

set(f4,'PaperPositionMode','auto');
set(f4,'PaperOrientation','landscape');
set(f4,'Position',[0 0 1920 1080]);
print(f4,'-painters','-dpdf','-bestfit','-r600',[fig_dir 'Participant Coefficient vs Module z-Score all networks.pdf']);

set(f5,'PaperPositionMode','auto');
set(f5,'PaperOrientation','landscape');
set(f5,'Position',[0 0 1920 1080]);
print(f5,'-painters','-dpdf','-bestfit','-r600',[fig_dir 'Guimera PC vs Module z-Score all networks.pdf']);
% print(f6,'-painters','-dpdf','-bestfit','-r600',[fig_dir 'Communicability all networks.pdf']);

set(f7,'PaperPositionMode','auto');
set(f7,'PaperOrientation','landscape');
set(f7,'Position',[0 0 1920 1080]);
print(f7,'-painters','-dpdf','-bestfit','-r600',[fig_dir 'Betweenness Centrality all networks.pdf']);

%--------------------------------------------------------------------------


function [Explore,randomOrdered100] = load_data()

%% Load Data
%Ask to load Zdenka or Adrian:
            waitfor(msgbox('Select the Explore saved data'));
            [FileName,PathName] = uigetfile('*.mat','Select the Explore saved data');
            f=fullfile(PathName,FileName);
            load(f);       
            wattsPath='C:\Users\aloe8475\Documents\PhD\GitHub\CODE\Analysis\Network Explore Functions';
            load([wattsPath 'Random_Ordered_Graphs_ALL_networks.mat']);
end

