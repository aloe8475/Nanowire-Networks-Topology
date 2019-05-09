function [Elout,AllFlag]=DetectThresoldCurrent(Sim,Electrodes,SimSets)
%% Detects if current thorugh any electrode exceeds thresold current
% and if so closes that electrode (by giving its openflag value 1)
MaxI=SimSets.Settings.MaxI;
AllFlag=0;
Elout=Electrodes;
if isequal(MaxI,0)
    return;
end




% 1) Take source names
it=height(Sim);
LastS=Sim(end,:);
ElList=LastS.Properties.VariableNames(contains(Sim.Properties.VariableNames,'ISource'));
TEl=length(ElList);
% 2) Compare currents with last I
for i=1:TEl
    ElCurr=LastS.(ElList{i});
    
    if  abs(ElCurr)>MaxI
        ixflag=find(contains(Electrodes.Name,ElList{i}(2:end)));
        of=Electrodes.OpenFlag{ixflag};
        of(it:end)=1;
        Elout.OpenFlag{ixflag}=of;
    end
    
end


ElFlag=cellfun(@(c) c(it),Elout.OpenFlag(contains(Electrodes.Name,'Source')));
if isequal(sum(ElFlag),TEl)
    AllFlag=1;
end




end