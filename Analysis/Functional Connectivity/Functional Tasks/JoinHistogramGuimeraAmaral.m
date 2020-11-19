load('PC_MZ_10V_5V_3V_2V.mat')

NLT10 = jh_cap([PC10{:}],[MZ10{:}],100,0,1,-5,5);
NLT5 = jh_cap([PC5{:}],[MZ5{:}],100,0,1,-5,5);
NLT3 = jh_cap([PC3{:}],[MZ3{:}],100,0,1,-5,5);
NLT2 = jh_cap([PC2{:}],[MZ2{:}],100,0,1,-5,5);

NLT2new=NLT2;
NLT3new=NLT3;
NLT5new=NLT5;
NLT10new=NLT10;

NLT2new(NLT2==0)=NaN;
NLT3new(NLT3==0)=NaN;
NLT5new(NLT5==0)=NaN;
NLT10new(NLT10==0)=NaN;
%%
imagescnan(NLT10new)
figure;
imagescnan(NLT2new)

% NLTMax = jh_cap(NLTPCMax,NLTMZMax,100,0,1,-5,5);
% NLTMin = jh_cap(NLTPCMin,NLTMZMin,100,0,1,-5,5);
%%
vmin = min([min(NLT10-NLT5), min(NLT10-NLT3), min(NLT5-NLT3),min(NLT10-NLT2),min(NLT5-NLT2),min(NLT3-NLT2)]);
vmax= max([max(NLT10-NLT5), max(NLT10-NLT3), max(NLT5-NLT3), max(NLT10-NLT2),max(NLT5-NLT2),max(NLT3-NLT2)]);
subplot(2,2,1);
a=NLT5-NLT10;
a(a==0)=NaN;
imagescnan(a)
title('5V - 10V');
caxis manual
caxis([vmin vmax]);
colorbar
%
subplot(2,2,2);
b=NLT3-NLT10;
b(b==0)=NaN;
imagescnan(b)
title('3V - 10V');
caxis manual
caxis([vmin vmax]);
colorbar
%
subplot(2,2,3);
c=NLT2-NLT10;
c(c==0)=NaN;
imagescnan(c)
title('2V - 10V');
caxis manual
caxis([vmin vmax]);
colorbar

%%
% %Loop through and get for each thing:
% load('Delta_MeanAcc_VSweep_WS.mat')
% 
% %BetaSweep
% for nn = 1:length(PC)
% jh_net_beta(:,:,nn) = jh_cap(PC(nn,:),MZ(nn,:),100,0,1,-6,6);
% end
% 
% for ii = 1:100
% for jj = 1:100
% [jh_mean_acc_beta(ii,jj),pval_mean_acc_beta(ii,jj)] = corr(squeeze(jh_net_beta(ii,jj,:)),transpose(MeanAcc));
% end
% end 
% for ii = 1:100
% for jj = 1:100
% [jh_delta_beta(ii,jj),pval_delta_beta(ii,jj)] = corr(squeeze(jh_net_beta(ii,jj,:)),transpose(Delta));
% end
% end 

%ASNs
clearvars PC MZ MeanAcc Delta jh_net jh_delta jh_mean_acc data PCtemp MZtemp MeanAcctemp Deltatemp
%Voltage Sweep:
networkType='NWN';
analType='Mod';%Mod 
if strcmp(networkType,'BA')
    onAmp={'0p2','0p5','0p75','1','1p25','1p5','2','5','10'};
    for i = 1:length(onAmp)
        data{i}=load(['Delta_MeanAcc_BA ' onAmp{i} 'V.mat']);
    end 
elseif strcmp(networkType,'NWN')
    onAmp={'0p2','0p5','0p75','1','1p25','1p5','2','5','10'};
    for i = 1:length(onAmp)
        if strcmp(analType, 'Mod')
            data{i}=load(['Delta_MeanAcc_NWN_9modules_' onAmp{i} 'V.mat']);
        else
            data{i}=load(['Delta_MeanAcc_NWN_AllDensities_' onAmp{i} 'V.mat']);
        end 
    end 
elseif strcmp(networkType,'HNW')
    if strcmp(analType, 'Mod')
        onAmp={'0p5','1','2'};
    else
        onAmp={'02','05','0p75','1','1p25','1p5','2','5','10'};
    end 
    for i = 1:length(onAmp)
        if strcmp(analType, 'Mod')
            data{i}=load(['Delta_MeanAcc_HNW Mod0p75 ' onAmp{i} 'V.mat']);
        else
            data{i}=load(['Delta_MeanAcc_V' onAmp{i} '.mat']);
        end 
    end 
end 

%HNW
for i = 1:length(data) %for each voltage 
    count=1;
    for j = 1:size(data{i}.MZ,1) %for each modularity
        for k = 1:size(data{i}.MZ,2) %for each network
            if ~all(data{i}.MZ{j,k}==0) %remove all networks with only MZ = 0
                PCtemp{i,j,k}=data{i}.PC{j,k};
                MZtemp{i,j,k}=data{i}.MZ{j,k};
                MeanAcctemp{i,j,k}=data{i}.MeanAcc(count);
                Deltatemp{i,j,k}=data{i}.Delta(count);
                count=count+1;
            end 
        end 
    end
end 
progressbar(0,0)

for i = 1:length(data) %each voltage
    progressbar([],0)
    for j = 1:size(data{i}.MZ,1) %for each modularity/density
        PC={PCtemp{i,j,:}};
        MZ={MZtemp{i,j,:}};
        MeanAcc=[MeanAcctemp{i,j,:}];
        Delta=[Deltatemp{i,j,:}];
        for nn = 1:length(PC)
            jh_net(:,:,nn) = jh_cap(PC{:,nn},MZ{:,nn},100,0,1,-6,10);
        end
        for ii = 1:100
            for jj = 1:100
                [jh_mean_acc{i,j}(ii,jj),pval_mean_acc{i,j}(ii,jj)] = corr(squeeze(jh_net(ii,jj,:)),transpose(MeanAcc));
            end
        end 
        for ii = 1:100
            for jj = 1:100
                [jh_delta{i,j}(ii,jj),pval_delta{i,j}(ii,jj)] = corr(squeeze(jh_net(ii,jj,:)),transpose(Delta));
            end
        end 
        progressbar([],j/20)

    end
    progressbar(i/9)

end 
 
if strcmp(networkType,'BA')
    BAdelta=jh_delta;
    BAma=jh_mean_acc;
elseif strcmp(networkType,'NWN')
    NWNdelta=jh_delta;
    NWNma=jh_mean_acc;
elseif strcmp(networkType,'HNW')
    HNWdelta=jh_delta;
    HNWma=jh_mean_acc;
end 

% %BA
% clearvars PC MZ MeanAcc Delta
% load('Delta_MeanAcc_VSweep_BA.mat')
% % load('Delta_MeanAcc_LessThanTwoVolts.mat')
% for nn = 1:length(PC)
% jh_net_BA(:,:,nn) = jh_cap(PC{:,nn},MZ{:,nn},100,0,1,-6,6);
% end
% 
% for ii = 1:100
% for jj = 1:100
% [jh_mean_acc_BA(ii,jj),pval_mean_acc_BA(ii,jj)] = corr(squeeze(jh_net_BA(ii,jj,:)),transpose(MeanAcc));
% end
% end 
% for ii = 1:100
% for jj = 1:100
% [jh_delta_BA(ii,jj),pval_delta_BA(ii,jj)] = corr(squeeze(jh_net_BA(ii,jj,:)),transpose(Delta));
% end
% end 
%%
%Plot betasweep
% for i =1:length(data)
%     a{i} = min(min([min(HNWdelta{i}),min(NWNdelta{i}),min(BAdelta{i})]),min([min(HNWma{i}),min(NWNma{i}),min(BAma{i})]));
%     b{i} = max(max([max(HNWdelta{i}),max(NWNdelta{i}),max(BAdelta{i})]),max([max(HNWma{i}),max(NWNma{i}),max(BAma{i})]));
% end
% bottom=min(min(min([NWNdelta{:}],[NWNma{:}])));
% top=max(max(max([NWNdelta{:}],[NWNma{:}])));
%% 
for i=1:size(NWNdelta,2) %for each density/modularity
    for j = 1:length(data) %for each voltage
        %plot nwn
        if j == 1
            f=figure('Renderer', 'painters', 'Position', [10 10 1400 1000]);
        end 
        set(0,'CurrentFigure',f)
        subplot(3,3,j)
        imagescnan(jh_delta{j,i})
        xlabel('PC');
        ylabel('MZ');
        caxis manual
        caxis([-1 1]);
        cb1=colorbar;
        ylabel(cb1,'r')
        % set(cb1,'YTick',[-250:50:50])
        title(['Delta (Standardized): V = ' onAmp{j}])
        if j == 1
            f2=figure('Renderer', 'painters', 'Position', [10 10 1400 1000]);
        end
        set(0,'CurrentFigure',f2)
        subplot(3,3,j)
        imagescnan(jh_mean_acc{j,i})
        title(['Mean Accuracy (Standardized): V = ' onAmp{j}])
        xlabel('PC');
        ylabel('MZ');
        caxis manual
        caxis([-1 1]);
        cb2=colorbar;
        ylabel(cb2,'r')
        % set(cb2,'YTick',[-250:50:50])
        if j == length(data)
            if strcmp(analType,'Mod')
                mod=[1,2,3,4,5,6];
                saveas(f,['../../../Data/Figures/Functional Connectivity/Node Level Comparisons/NWN/Voltage Sweep/Modularity Comparison/' networkType ' VSweep Delta (NLT-MC) PC+MZ Correlations mod_' int2str(mod(i)) '.png']);
                saveas(f2,['../../../Data/Figures/Functional Connectivity/Node Level Comparisons//NWN/Voltage Sweep/Modularity Comparison/' networkType ' VSweep Mean Accuracy (NLT+MC) PC+MZ Correlations mod_' int2str(mod(i)) '.png']);
            else
                avgDeg=[4,   5,   5,   6,   7,   8,   9,  10,  12,  14,  17, 20,  27,  34,  45,  68,  98, 170, 238, 285];
                saveas(f,['../../../Data/Figures/Functional Connectivity/Node Level Comparisons/' networkType ' Networks/Voltage Sweep/Density Comparison/' networkType ' VSweep Delta (NLT-MC) PC+MZ Corr avgDeg_' int2str(avgDeg(i)) '.png']);
                saveas(f2,['../../../Data/Figures/Functional Connectivity/Node Level Comparisons/' networkType ' Networks/Voltage Sweep/Density Comparison/' networkType ' VSweep Mean Accuracy (NLT+MC) PC+MZ Corr avgDeg_' int2str(avgDeg(i)) '.png']);
                close all
            end
        end
    %         saveas(f2,['../../../Data/Figures/Functional Connectivity/Node Level Comparisons/NWN V' onAmp{i} ' Mean Accuracy (NLT+MC) PC+MZ Correlations.png']);
    % %plot BA
    % f=figure();
    % imagescnan(jh_delta_BA)
    % xlabel('PC');
    % ylabel('MZ');
    % caxis manual
    % caxis([bottom top]);
    % cb1=colorbar;
    % ylabel(cb1,'r')
    % % set(cb1,'YTick',[-250:50:50])
    % title('Delta (NLT - MC) Standardized Accuracy: BA Correlations')
    % f2=figure();
    % imagescnan(jh_mean_acc_BA)
    % title('Mean Accuracy NLT + MC: BA Correlations')
    % xlabel('PC');
    % ylabel('MZ');
    % caxis manual
    % caxis([bottom top]);
    % cb2=colorbar;
    % ylabel(cb2,'r')
    % % set(cb2,'YTick',[-250:50:50])
    % saveas(f,'../../../Data/Figures/Functional Connectivity/Node Level Comparisons/BA Delta PCOEFF + MZ Correlations - 2V.png');
    % saveas(f2,'../../../Data/Figures/Functional Connectivity/Node Level Comparisons/BA Mean Accuracy PCOEFF + MZ Correlations - 2V.png');
    end 
end 
%%
pval_delta_beta_threshold05=double(pval_delta_beta<0.05 & pval_delta_beta>0.001);
pval_delta_beta_threshold001=double(pval_delta_beta<0.001);
pval_delta_ba_threshold05=double(pval_delta_BA<0.05 & pval_delta_BA>0.001);
pval_delta_ba_threshold001=double(pval_delta_BA<0.001);
pval_delta_threshold05=double(pval_delta<0.05 & pval_delta>0.001);
pval_delta_threshold001=double(pval_delta<0.001);

pval_delta_beta_threshold001(pval_delta_beta_threshold001==1)=2;
pval_delta_beta_threshold001(pval_delta_beta_threshold05==1)=1;
pval_delta_ba_threshold001(pval_delta_ba_threshold001==1)=2;
pval_delta_ba_threshold001(pval_delta_ba_threshold05==1)=1;
pval_delta_threshold001(pval_delta_threshold001==1)=2;
pval_delta_threshold001(pval_delta_threshold05==1)=1;


pval_mean_acc_beta_threshold05=double(pval_mean_acc_beta<0.05 & pval_mean_acc_beta>0.001);
pval_mean_acc_beta_threshold001=double(pval_mean_acc_beta<0.001);
pval_mean_acc_threshold05=double(pval_mean_acc<0.05 & pval_mean_acc>0.001);
pval_mean_acc_threshold001=double(pval_mean_acc<0.001);
pval_mean_acc_ba_threshold05=double(pval_mean_acc_BA<0.05 & pval_mean_acc_BA>0.001);
pval_mean_acc_ba_threshold001=double(pval_mean_acc_BA<0.001);

pval_mean_acc_beta_threshold001(pval_mean_acc_beta_threshold001==1)=2;
pval_mean_acc_beta_threshold001(pval_mean_acc_beta_threshold05==1)=1;
pval_mean_acc_ba_threshold001(pval_mean_acc_ba_threshold001==1)=2;
pval_mean_acc_ba_threshold001(pval_mean_acc_ba_threshold05==1)=1;
pval_mean_acc_threshold001(pval_mean_acc_threshold001==1)=2;
pval_mean_acc_threshold001(pval_mean_acc_threshold05==1)=1;


%%
%Plot betasweep pvals
f=figure();
bottom = 0; %min(min(min(min(pval_delta_beta),min(pval_delta)),min(min(pval_mean_acc_beta),min(pval_mean_acc))));
top  = 2;%max(max(max(max(pval_delta_beta),max(pval_delta)),max(max(pval_mean_acc_beta),max(pval_mean_acc))));
h=imagescnan(pval_delta_beta_threshold001);
imlegend(h,[1,2],{'p<0.05','p<0.001'})
xlabel('PC');
ylabel('MZ');
caxis manual
caxis([bottom top]);
% colormap(f,'prism')
% cb1=colorbar;
% ylabel(cb1,'r')
% set(cb1,'Ticks',[1,2])
% set(cb1,'TickLabels',[0.01,0.05])
title('Delta Watts-Strogatz Correlation PVal')
f2=figure();
h2=imagescnan(pval_mean_acc_beta_threshold001);
imlegend(h2,[1,2],{'p<0.05','p<0.001'})

title('Mean Accuracy Watts-Strogatz Correlation PVal')
xlabel('PC');
ylabel('MZ');
caxis manual
caxis([bottom top]);

% colormap(f2,'prism')
% cb2=colorbar;
% ylabel(cb2,'r')
% set(cb2,'Ticks',[1,2])
% set(cb2,'TickLabels',[0,0.05])
saveas(f,'../../../Data/Figures/Functional Connectivity/Node Level Comparisons/WS VSweep Delta (NLT-MC) PC+MZ Correlations PVals.png');
saveas(f2,'../../../Data/Figures/Functional Connectivity/Node Level Comparisons/WS VSweep Mean Accuracy (NLT+MC) PC+MZ Correlations PVals.png');

%plot BA
f=figure();
h=imagescnan(pval_delta_ba_threshold001)
imlegend(h,[1,2],{'p<0.05','p<0.001'})
xlabel('PC');
ylabel('MZ');
caxis manual
caxis([bottom top]);
% colormap(f,'prism')
% cb1=colorbar;
% ylabel(cb1,'r')
% set(cb1,'Ticks',[1,2])
% set(cb1,'TickLabels',[0,0.05])
title('Delta BA Correlation PVal')
f2=figure();
h2=imagescnan(pval_mean_acc_ba_threshold001)
imlegend(h2,[1,2],{'p<0.05','p<0.001'})
title('Mean Accuracy BA Correlation PVal')
xlabel('PC');
ylabel('MZ');
caxis manual
caxis([bottom top]);
% colormap(f2,'prism')
% cb2=colorbar;
% ylabel(cb2,'r')
% set(cb2,'Ticks',[1,2])
% set(cb2,'TickLabels',[0,0.05])
saveas(f,'../../../Data/Figures/Functional Connectivity/Node Level Comparisons/BA VSweep Delta (NLT-MC) PC+MZ Correlations PVals.png');
saveas(f2,'../../../Data/Figures/Functional Connectivity/Node Level Comparisons/BA VSweep Mean Accuracy (NLT+MC) PC+MZ Correlations PVals.png');


%plot nwn
f=figure();
h=imagescnan(pval_delta_threshold001)
imlegend(h,[1,2],{'p<0.05','p<0.001'})
xlabel('PC');
ylabel('MZ');
caxis manual
caxis([bottom top]);
% colormap(f,'prism')
% cb1=colorbar;
% ylabel(cb1,'r')
% set(cb1,'Ticks',[1,2])
% set(cb1,'TickLabels',[0,0.05])
title('Delta NWN Correlation PVal')
f2=figure();
h2=imagescnan(pval_mean_acc_threshold001)
imlegend(h2,[1,2],{'p<0.05','p<0.001'})
title('Mean Accuracy NWN Correlation PVal')
xlabel('PC');
ylabel('MZ');
caxis manual
caxis([bottom top]);
% colormap(f2,'prism')
% cb2=colorbar;
% ylabel(cb2,'r')
% set(cb2,'Ticks',[1,2])
% set(cb2,'TickLabels',[0,0.05])
saveas(f,'../../../Data/Figures/Functional Connectivity/Node Level Comparisons/NWN VSweep Delta (NLT-MC) PC+MZ Correlations PVals.png');
saveas(f2,'../../../Data/Figures/Functional Connectivity/Node Level Comparisons/NWN VSweep Mean Accuracy (NLT+MC) PC+MZ Correlations PVals.png');

