function [mSour,mOff]=ControlAverager(currTimeO,currSour,currOff,SelSour,SelOff,AveTime,sliderT)



TimeIndex=round(sliderT.Value);
T=sliderT.UserData(TimeIndex);
TimeIndex=find(currTimeO<=T,1);

if isempty(TimeIndex)
    TimeIndex=1;
end

mSour=zeros(1,9);
mOff=zeros(1,9);
switch AveTime
    case 'UptoR'
        mSour(SelSour)=mean(currSour(1:TimeIndex,SelSour),1);
        mOff(SelOff)=mean(currOff(1:TimeIndex,SelOff),1);
    case 'NoAveR'
        mSour(SelSour)=currSour(TimeIndex,SelSour);
        mOff(SelOff)=currOff(TimeIndex,SelOff);
    case 'AveAllR'
        mSour(SelSour)=mean(currSour(:,SelSour),1);
        mOff(SelOff)=mean(currOff(:,SelOff),1);       
end
end