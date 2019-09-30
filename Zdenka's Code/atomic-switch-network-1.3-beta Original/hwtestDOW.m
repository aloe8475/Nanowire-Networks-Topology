global t;
global a;
global plotpsd;
global dutyratio;
global ap;
global Nodes
for Nodes = 1:10
contactn=randperm(100,3);
 
% load('handwritedata.mat');load('handwirteresult.mat');
Data=dlmread('C:\Users\aloe8475\Documents\PhD\GitHub\CODE\Analysis\Classification\Time Serieis\Dow_Jones_Percentage_Change1.csv');

% hwresult=Data;
% % global b;
% % global c;
% % if 1
% %     recurrent=[];
% N=hwresult(1:end);
% % M=(hwdata(1:6000));
% result=zeros(length(N),10);
% for i=1:length(N)
% r=N(i);result(i,r+1)=1;
% end
% for ap=0.5:0.4:2
% for b=1:length(M)
    [testcurrent testconduct]=MonolithicDemo(contactn);
    testcurrent(end)
    recurrent(Nodes,:)=testcurrent;
    reconduct(Nodes,:)=testconduct;
end
%     recurent(b,:)=Output.networkCurrent(c,:);
% end
% recurent(:,end)=[];
save(['DOW_current_conduct1.mat'],'recurrent','reconduct')
% save(['DOW_result' num2str(ap) '.mat'],'result')
% end
% end 
% if 0
% N=hwresult(501:1000);
% M=hwdata(501:1000);
% Ntest=zeros(length(N),1);Ntarget=Ntest;
% result=zeros(length(N),10);
%  for i=1:length(N)
% r=N(i);result(i,r+1)=1;
% end
% for b=1:length(M)
%     MonolithicDemokevin;
%     recurent(b,:)=Output.networkCurrent(end,:);
% %     re=recurent(b,:);re(:,c)=[];
%     cc=ann28(testcurrent);
% [maxVal maxInd] =max(cc);
% Ntest(b)=maxInd-1;
% for i=1:length(M)
% [maxVal maxInd]=max(result(i,:));
% Ntarget(i)=maxInd-1;
% end
% end
% accuracy=sum(Ntest==Ntarget)/length(Ntest);
% end
