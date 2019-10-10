
addpath(genpath(pwd));
contactn=[randperm(100,29)];
load('handwirteresult.mat');load('hwdata.mat');
hwresult=c;
if 1
    testcurrent=[];
    recurrent=[];
N=hwresult(1:60000);
N=reshape(N,[10,6000]);
c=reshape(c,[10,6000]);

result=zeros(length(N),10);
for i=1:length(N)
r=N(i);result(i,r+1)=1;
end

numworkers=6;

pool=parpool(numworkers);
for ap=0.6:0.2:1.6
parfor w=1:numworkers
for    b=1:size(M,2)
[testcurrent]=MonolithicDemoMNIST(contactn,b,c(w,:),M(w,:),N(w,:),ap);
    recurrent{w}=testcurrent;
	end
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
