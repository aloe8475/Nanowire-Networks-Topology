function ref=OpenTFAssoc(fpath)
[Out,~]=TDMS_readTDMSFile(fpath);
%%

%Extract relevant information from tdms
%channel info

[textdata,~]=prepareMetaData(Out.propNames{1},Out.propValues{1});



Out.data(cellfun('isempty',Out.data))=[];
data=Out.data;

indexoutput=Out.chanIndices{1};
indexcurrent=Out.chanIndices{2};




 % Get data current or voltage 
length_Data=length(Out.data{indexcurrent(1)});
curr=zeros(length_Data,length(indexcurrent));
j=1;
for i=(indexcurrent-2)'
    curr(:,j)=data{i};
    j=j+1;
end
channels=Out.chanNames{2}';
   %subdivide in cells
 list999=find(curr(:,1)==999);
 currdata=cell(1,length(list999));
 timedata=cell(1,length(list999));
 Rate=str2double(textdata.DefFreqAcq);
 for i=1:(length(list999)-1)
     currdata{i}=curr(list999(i)+1:list999(i+1)-1,:);
     timedata{i}=linspace(0,(1/Rate)*(length(currdata{i}(:,1))),length(currdata{i}(:,1)));
 end
 currdata{i+1}=curr(list999(i+1)+1:end,:);
 timedata{i+1}=linspace(0,(1/Rate)*(length(currdata{i+1}(:,1))),length(currdata{i+1}(:,1)));

 
 %% GET OUTPUT INFO

length_Out=length(Out.data{indexoutput(2)});
curro=zeros(length_Out,length(indexoutput));
j=1; 
for i=(indexoutput-3)'
    curro(:,j)=data{i};
    j=j+1;
end

%for the moment I'll just assume 9 source and 9 offset
channelsSource={'Source0','Source1','Source2','Source3','Source4','Source5','Source6','Source7','Source8'}';
channelsOff={'Offset0','Offset1','Offset2','Offset3','Offset4','Offset5','Offset6','Offset7','Offset8'}';
   %subdivide in cells
 listNaN=find(isnan(curro(:,2)));
 currsource=cell(1,length(listNaN));
 curroff=cell(1,length(listNaN));;
 timeout=cell(1,length(listNaN));

 for i=1:(length(listNaN)-1)
     currsource{i}=curro(listNaN(i)+1:listNaN(i+1)-1,2:10);
     curroff{i}=curro(listNaN(i)+1:listNaN(i+1)-1,11:19);
     timeout{i}=curro(listNaN(i)+1:listNaN(i+1)-1,1)-curro(listNaN(i),1);
 end
 currsource{i+1}=curro(listNaN(i+1)+1:end,2:10);
 curroff{i+1}=curro(listNaN(i)+1:listNaN(i+1)-1,11:19);
 timeout{i+1}=curro(listNaN(i+1)+1:end,1)-curro(listNaN(i+1),1);

 
% Get Additional Info
ControlInfo=textdata.CompleteTrainTestInformation;
listTrTs=find(ControlInfo=='T');
PatternInfo=cell(1,length(listTrTs));
DIDOlist=cell(1,length(listTrTs));
TrTs=cell(1,length(listTrTs));
DinNumb=zeros(1,length(listTrTs));
DoutNumb=zeros(1,length(listTrTs));
DtarNumb=zeros(1,length(listTrTs));
for i=1:length(listTrTs)
    aux.StrPat=ControlInfo(listTrTs(i):listTrTs(i)+58);
    [dido.digin,dido.digout,dido.target]=split_patstr(aux.StrPat);
    TrTs{i}=aux.StrPat(1:2);
    DinNumb(i)=bin2dec(fliplr(dido.digin));
    DoutNumb(i)=bin2dec(fliplr(dido.digout));
    DtarNumb(i)=bin2dec(fliplr(dido.target));
    PatternInfo{i}=aux.StrPat;
    DIDOlist{i}=dido;    
end


ref.textdata=textdata;
ref.Filename=textdata.name;

ref.channels=channels';
ref.timedata=timedata;
ref.currdata=currdata;

ref.channelsSource=channelsSource;
ref.channelsOff=channelsOff;
ref.currsource=currsource;
ref.curroff=curroff;
ref.timeout=timeout;


ref.PatternInfo=PatternInfo;
ref.DIDOlist=DIDOlist;
ref.DinNumb=DinNumb;
ref.DoutNumb=DoutNumb;
ref.DtarNumb=DtarNumb;
ref.TrTs=TrTs;





end