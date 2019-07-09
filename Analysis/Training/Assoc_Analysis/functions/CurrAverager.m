function [mData]=CurrAverager(~,currData,SelChan,AveTime,sliderT)



TimeIndex=round(sliderT.Value);


mData=zeros(1,9);

switch AveTime
    case 'UptoR'
        mData(SelChan)=mean(currData(1:TimeIndex,SelChan),1);
    case 'NoAveR'
        mData(SelChan)=currData(TimeIndex,SelChan);
    case 'AveAllR'
        mData(SelChan)=mean(currData(:,SelChan),1);
        
end
end