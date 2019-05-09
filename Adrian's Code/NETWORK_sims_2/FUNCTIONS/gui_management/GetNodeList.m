function [NodeList,UserData]=GetNodeList(Sim)
   Layout=Sim.SelLayout;
   [i,j]=find(Layout.AdjMat);
   NodeList=cell(length(i),1);
   for k=1:length(i)
       NodeList{k}=strcat('Node_',num2str(j(k)),'_',num2str(i(k)));
   end
   Electrodes=Sim.Settings.ElectrodesInfo;
   i=length(Electrodes);
   ElList=cell(i,1);
   UserData.SelectedNodes=i;
   for k=1:i
       ElList{k}=strcat('I',Electrodes{k}.Name);
   end
   NodeList=vertcat(ElList,NodeList);
end