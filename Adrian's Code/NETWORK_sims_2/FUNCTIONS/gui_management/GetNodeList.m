function [NodeList,UserData]=GetNodeList(Sim,varargin)
Layout=Sim.SelLayout;
[i,j]=find(Layout.AdjMat);
NodeList=cell(length(i),1);
for k=1:length(i)
    NodeList{k}=strcat('Node_',num2str(j(k)),'_',num2str(i(k)));
end
if ~isempty(varargin) & varargin == 'z' 
    Electrodes=[Sim.Electrodes.PosIndex]; %NEED TO FIX FOR ZDENKA to ADRIAN CODE
    i=length(Electrodes);
    ElList=cell(i,1);
    UserData.SelectedNodes=i;
    Names={Sim.Electrodes.Name};
    for k=1:i
        ElList{k}=strcat('I',Names{k});
    end
NodeList=vertcat(ElList,NodeList);

else
    Electrodes=[Sim.Electrodes.PosIndex]; %NEED TO FIX FOR ZDENKA to ADRIAN CODE
    i=length(Electrodes);
    ElList=cell(i,1);
    UserData.SelectedNodes=i;
    Names={Sim.Electrodes.Name};
    for k=1:i
        ElList{k}=strcat('I',Names{1}{k});
    end
NodeList=vertcat(ElList,NodeList);
end