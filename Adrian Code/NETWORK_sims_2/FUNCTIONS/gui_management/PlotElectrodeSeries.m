function PlotElectrodeSeries(currAx,NodeList,Sim,type,status,timeInd,varargin)
%% plot time series of network


els=cellstr(Sim.Electrodes.Name(contains(Sim.Electrodes.Name,'Source')));
eld=cellstr(Sim.Electrodes.Name(contains(Sim.Electrodes.Name,'Drain')));

Ns=ceil(sqrt(length(els)));
Nd=ceil(sqrt(length(eld)));
cla(currAx);

if contains(type,'Cond')
    stype='Conductance';
else
    stype='Current';
end

if contains(type,'S')
    ElSel=els;
    N=Ns;
else
    ElSel=eld;
    N=Nd;
end
ElSel=strcat('I',ElSel);

[~,~,~,cy,~]=ExtractCurrentDataFromSim(ElSel,Sim.Data,stype,timeInd);


if length(cy)>3
    numel=N^2;
    if ~isequal(length(cy),numel)
        cy=[cy NaN(1,numel-length(cy))];
    end
    Cy=reshape(cy,[N,N]);
else
    Cy=cy;
end


image(currAx,abs(Cy),'CDataMapping','scaled');
colormap(currAx,gcurrmap);
colorbar(currAx);
caxis(currAx,[Sim.SimInfo.MinI Sim.SimInfo.MaxI]);
im=currAx.Children;
X=im.XData;
Y=im.YData;
if length(X)>1
    wx=(X(2)-X(1))/(size(im.CData,2)-1);
else
    wx=X;

end
if length(Y)>1
    wy=(Y(2)-Y(1))/(size(im.CData,2)-1);
else
    wy=Y;
end
x=[];y=[];
for i=1:X(end)
    for j=1:Y(end)
        x=[x i*wx];
        y=[y j*wy];
    end
end

if length(x)>length(cellstr(ElSel))
    dif=length(x)-length(cellstr(ElSel));
    of=length(cellstr(ElSel));
    for k=of:of+dif
        ElSel{k}='-';
    end
end

text(currAx,x,y,ElSel,'FontSize',10,'FontWeight','bold','HorizontalAlignment','center','Color','yellow');












end
