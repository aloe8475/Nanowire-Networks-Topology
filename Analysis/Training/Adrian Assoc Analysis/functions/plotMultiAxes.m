function h=plotMultiAxes(handles,Sel,SelChan,SelSource)
reflist=handles.refList{1};
Sel=Sel(1);
currTime=reflist.timedata{Sel};
currData=reflist.currdata{Sel};
currTimeO=reflist.timeout{Sel};
currSour=reflist.currsource{Sel};
currOff=reflist.curroff{Sel};
SelSour=SelSource(SelSource<10);
SelOff=SelSource(SelSource>=10);SelOff=SelOff-9;

[mData]=CurrAverager(currTime,currData,SelChan,...
    handles.AveTime,handles.TimeSlider);
[mSour,mOff]=ControlAverager(currTimeO,currSour,currOff,SelSour,SelOff,...
    handles.AveTime,handles.TimeSlider);


currData=currData(:,SelChan);

if ~isempty(SelSour)
    currSour=currSour(:,SelSour);
    rm=repmat([1 0 0],length(SelSour),1);
else
    currSour=[];rm=[];
end
if ~isempty(SelOff)
    currOff=currOff(:,SelOff);
    bm=repmat([0 0 1],length(SelOff),1);
else
    currOff=[];bm=[];
end
currDOut=[currSour,currOff];
corder=[rm;bm];







hold(handles.axesT,'off');
if ~isempty(handles.timeLineaxT)
    X=handles.timeLineaxT.XData;
    Y=handles.timeLineaxT.YData;
    plot(handles.axesT,X,Y,'Color','r','Linewidth',0.9);
    handles.timeLineaxT=handles.axesT.Children;
    hold(handles.axesT,'on');
end
plot(handles.axesT,currTime,currData);
mT=max(currTime);
set(handles.axesT,'XLim',[0 mT]);


hold(handles.axesControl,'off');
if ~isempty(handles.timeLineaxCtrl)
    X=handles.timeLineaxCtrl.XData;
    Y=handles.timeLineaxCtrl.YData;
    plot(handles.axesControl,X,Y,'Color','r','Linewidth',0.9);
    handles.timeLineaxCtrl=handles.axesControl.Children;
    hold(handles.axesControl,'on');
end

for i=1:size(currDOut,2)
    plot(handles.axesControl,currTimeO,currDOut(:,i),'Color',corder(i,:));
    hold(handles.axesControl,'on');
end

mT=max(currTimeO);
set(handles.axesControl,'XLim',[0 mT]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
colormap(handles.axesCurrent,handles.defColormap);
ColorMapCurrent=reshape(mData,3,3)';
imagesc(handles.axesCurrent,ColorMapCurrent,'CDataMapping',handles.Scaling);
if handles.Scaling=='scaled'
    caxis(handles.axesCurrent,[-10 10]);
end
colorbar(handles.axesCurrent,'Location','westoutside');
set(handles.axesCurrent,'XTick',[]);
set(handles.axesCurrent,'YTick',[]);


colormap(handles.axesOff,handles.defColormap);
ColorMapCurrent=reshape(mOff,3,3)';
imagesc(handles.axesOff,ColorMapCurrent,'CDataMapping',handles.Scaling);
if handles.Scaling=='scaled'
    caxis(handles.axesOff,[-10 10]);
end
caxis(handles.axesOff,[-10 10]);
set(handles.axesOff,'XTick',[]);
set(handles.axesOff,'YTick',[]);


colormap(handles.axesSource,handles.defColormap);
ColorMapCurrent=reshape(mSour,3,3)';
imagesc(handles.axesSource,ColorMapCurrent,'CDataMapping',handles.Scaling);
if handles.Scaling=='scaled'
    caxis(handles.axesSource,[-10 10]);
end
set(handles.axesSource,'XTick',[]);
set(handles.axesSource,'YTick',[]);
colorbar(handles.axesSource,'Location','eastoutside');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
boolarr=ToBoolArray(reflist.DIDOlist{Sel}.digin);
ColorMapDig=reshape(boolarr,3,3)';
imagesc(handles.axesDigInput,ColorMapDig,'CDataMapping','scaled');
caxis(handles.axesDigInput,[0 1]);
set(handles.axesDigInput,'XTick',[]);
set(handles.axesDigInput,'YTick',[]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
boolarr=ToBoolArray(reflist.DIDOlist{Sel}.digout);
ColorMapDig=reshape(boolarr,3,3)';
imagesc(handles.axesDigOutput,ColorMapDig,'CDataMapping','scaled');
caxis(handles.axesDigOutput,[0 1]);
set(handles.axesDigOutput,'XTick',[]);
set(handles.axesDigOutput,'YTick',[]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
boolarr=ToBoolArray(reflist.DIDOlist{Sel}.target);
ColorMapDig=reshape(boolarr,3,3)';
imagesc(handles.axesTarget,ColorMapDig,'CDataMapping','scaled');
caxis(handles.axesTarget,[0 1]);
set(handles.axesTarget,'XTick',[]);
set(handles.axesTarget,'YTick',[]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=handles;



end