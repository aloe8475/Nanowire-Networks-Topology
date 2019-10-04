contactn=[randperm(100,29)];
load('handwritedata.mat');load('handwirteresult.mat');
hwdata=b;
hwresult=c;
if 1
    recurrent=[];
N=hwresult(1:60000);
M=(hwdata(1:60000));
result=zeros(length(N),10);
for i=1:length(N)
r=N(i);result(i,r+1)=1;
end
pool=parpool(10);
for ap=0.6:0.2:1.6
parfor b=1:length(M)
   [testcurrent]=MonolithicDemoMNIST(contactn,b,c,M,N,ap);
    recurrent=[recurrent;testcurrent];
%     recurent(b,:)=Output.networkCurrent(c,:);
end
% recurent(:,end)=[];
save(['60000hwdata_Amp' num2str(ap) '_700nw.mat'],'recurrent')
save(['60000hwdata_Result' num2str(ap) '_700nw.mat'],'result')
end
end 
if 0
N=hwresult(501:1000);
M=hwdata(501:1000);
Ntest=zeros(length(N),1);Ntarget=Ntest;
result=zeros(length(N),10);
 for i=1:length(N)
r=N(i);result(i,r+1)=1;
end
for b=1:length(M)
    MonolithicDemokevin;
    recurent(b,:)=Output.networkCurrent(end,:);
%     re=recurent(b,:);re(:,c)=[];
    cc=ann28(testcurrent);
[maxVal maxInd] =max(cc);
Ntest(b)=maxInd-1;
for i=1:length(M)
[maxVal maxInd]=max(result(i,:));
Ntarget(i)=maxInd-1;
end
end
accuracy=sum(Ntest==Ntarget)/length(Ntest);
end
