function Kdata=Kaverager(~,currData,SelChan,AveTime,sliderT)



TimeIndex=round(sliderT.Value);




switch AveTime
    case 'UptoR'
        Kdata=zeros(length(1:TimeIndex),9);
        Kdata(:,SelChan)=currData(1:TimeIndex,SelChan);
    case 'NoAveR'
        Kdata=zeros(1,9);
        Kdata(1,SelChan)=currData(TimeIndex,SelChan);
    case 'AveAllR'
        Kdata=zeros(size(currData));
        Kdata(:,SelChan)=currData(:,SelChan);
        
end
end
