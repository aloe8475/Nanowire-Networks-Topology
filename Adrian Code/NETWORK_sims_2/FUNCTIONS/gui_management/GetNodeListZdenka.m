function [NodeList,UserData]=GetNodeListZdenka(Sim,network_load)
Layout=Sim.SelLayout;
[i,j]=find(Layout.AdjMat);
NodeList=cell(length(i),1);
for k=1:length(i)
    NodeList{k}=strcat('Node_',num2str(j(k)),'_',num2str(i(k)));
end
if network_load=='a'
    Electrodes=Sim.Settings.ElectrodesInfo; %NEED TO FIX FOR ZDENKA to ADRIAN CODE
    i=length(Electrodes);
    ElList=cell(i,1);
    UserData.SelectedNodes=i;
    for k=1:i
        ElList{k}=strcat('I',Electrodes{k}.Name);
    end
else
    Electrodes={Sim.Electrodes.Name}; %NEED TO FIX FOR ZDENKA to ADRIAN CODE
    i=length(Electrodes);
    ElList=cell(i,1);
    UserData.SelectedNodes=i;
    for k=1:i
        ElList{k}=strcat('I',Electrodes{k});
    end
end
NodeList=vertcat(ElList,NodeList);
end