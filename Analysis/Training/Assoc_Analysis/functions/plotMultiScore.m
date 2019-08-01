function h=plotMultiScore(handles,Sel,SelChan)
reflist=handles.refList{1};
Dtar=reflist.DtarNumb;
ScoreAlg=handles.ScoringPop.String{handles.ScoringPop.Value};
Kchoose=handles.Kfitgroup.SelectedObject.Tag;
Knumb=str2double(handles.MaxKEdit.String);
if (Knumb<1)
    Knumb=1;
end

j=1;
mData=cell(1,length(Sel));
KData=cell(1,length(Sel));
for i=Sel
currData=reflist.currdata{i};
currTime=reflist.timedata{i};
mData{j}=CurrAverager(currTime,currData,SelChan,...
    handles.AveTime,handles.TimeSlider);
KData{j}=Kaverager(currTime,currData,SelChan,...
    handles.AveTime,handles.TimeSlider);
Tar(j)=Dtar(i);
j=j+1;
end

%If we want to save the currents
path='D:\alon_\Research\PhD\CODE\Data\Associative Training Data\Output Data\';
save([path reflist.Filename '.mat'],'reflist');

if isequal(ScoreAlg,'Kmeans')
    [seq,col,centroids]=scorer(KData,Tar,ScoreAlg,Kchoose,Knumb);
    ScoreList=cell(2,1);
    Cstr=RemoveWhiteSpaces(num2str(sort(centroids)),';');
    Tstr=RemoveWhiteSpaces(num2str(sort(Tar)),';');
    
    ScoreList{1}=strcat(ScoreAlg,'_centroids :',Cstr);
    ScoreList{2}=strcat(ScoreAlg,'_targets :',Tstr);
    
    scatter(handles.axesAux,seq,seq,'CData',col);
    boolarr=zeros(length(centroids),9);
    for i=1:length(centroids)
        boolarr(i,:)=NumbToBoolArray(round(centroids(i)));
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    imagesc(handles.axesCentroid,boolarr,'CDataMapping','scaled');
    caxis(handles.axesCentroid,[0 1]);
    set(handles.axesCentroid,'XTick',[]);
    set(handles.axesCentroid,'YTick',[]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
else
    [scores,~,~]=scorer(mData,Tar,ScoreAlg,[],[]);
    ScoreList=cell(length(scores),1);
    for i=1:length(scores)
        ScoreList{i}=strcat(ScoreAlg,':',num2str(scores(i)),'/ Tar:',num2str(Tar(i)));
    end
    scatter(handles.axesAux,categorical(Tar),scores);
    
    
end


handles.axesAux.Box='On';
handles.axesAux.Color=[0.5 0.5 0.5];
handles.ScoreText.String=ScoreList;














%%%%%%%%%MULTIMATRIX PLOTTING
%only first matrix of selected will be plot
mData=mData{1};
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





Sel=Sel(1);
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