function PlotNetworkLayout(currAx,NodeList,Sim,type,Status,IndexTime,varargin)

%%retrieve sim info
Layout=Sim.SelLayout;

x1=diag(Layout.X1);
x2=diag(Layout.X2);
y1=diag(Layout.Y1);
y2=diag(Layout.Y2);
X=full([x1' ; x2']);
Y=full([y1' ; y2']);
Gr=Layout.SelGraph;
[~,~,Cx]=find(Layout.CX);
[~,~,Cy]=find(Layout.CY);
Adj=triu(Layout.AdjMat);
NumEl=height(Sim.Electrodes);


IdxEl=Sim.Electrodes.PosIndex(NodeList.Value(NodeList.Value<=NumEl));

[SelEl,NodeSelMat,NodeSelList]=SelFromNodeList(length(Adj),NodeList);





%% network layout
if ~isequal(Status,'redraw') && ~isequal(type,'VoltageLayout')
    cla(currAx);
    PlotNetworkAux(currAx,X,Y,Cx,Cy,'all');
end
hold(currAx,'on');
%% electrodes are plot if selected in nodelist
if  ~isempty(SelEl) && ~isequal(type,'VoltageLayout')
    Xe=X(:,IdxEl);Ye=Y(:,IdxEl);
    Cxe=(x2(IdxEl)+x1(IdxEl))./2;Cye=(y1(IdxEl)+y2(IdxEl))./2;
    PlotNetworkAux(currAx,Xe,Ye,Cxe,Cye,'electrodes',SelEl,...
        NodeList.UserData.Annotate);
else
    %erase drawn electrodes if any
    sc=findobj(currAx,'LineWidth',2.2);
    tx=findobj(currAx,'Type','text');
    if ~isempty(sc) delete(sc); end %#ok<*SEPEX>
    if ~isempty(tx) delete(findobj(tx,'Color','green')); end
end

%% additional information to draw on the network
hold(currAx,'on');
switch type
    case 'SelEdges'  %%selected nodes
        if ~isequal(nnz(NodeSelMat),0)
            CxSel=full(Layout.CX(NodeSelMat~=0));
            CySel=full(Layout.CY(NodeSelMat~=0));
            PlotNetworkAux(currAx,X,Y,CxSel,CySel,'nodes',...
                NodeSelList,NodeList.UserData.Annotate);
        else
            %erase drawn nodes if any
            sc=findobj(currAx,'MarkerFaceColor','white');
            tx=findobj(currAx,'Type','text');
            if ~isempty(sc) delete(sc); end %#ok<*SEPEX>
            if ~isempty(tx) delete(findobj(tx,'Color','white')); end
        end
    case 'ShortestPath' % shortest path between every pair of selected electrodes
        pt=findobj(currAx,'Color','black');
        if ~isempty(pt) delete(pt); end
        if length(IdxEl)>=2
            combs=nchoosek(IdxEl,2); %combs between pair of electrodes
            path=[];
            for i=1:size(combs,1)
                path=[path shortestpath(Gr,combs(i,1),combs(i,2))]; %#ok<AGROW>
            end
            path=unique(path);
            X=X(:,path);Y=Y(:,path);
            PlotNetworkAux(currAx,X,Y,Cx,Cy,'path','k');
        end
    case 'VoltageLayout' %voltage map of network
        vlist=Sim.Data.Voltages{IndexTime};
        %interpolate voltage and color
        v=linspace(Sim.SimInfo.MinV,Sim.SimInfo.MaxV,10*length(vlist));
        cmap=hot(length(v));
        c=interp1(full(v),cmap,full(vlist));
        c(isnan(c))=0;
        PlotNetworkAux(currAx,X,Y,Cx,Cy,'volt',c);
        clim=[Sim.SimInfo.MinV Sim.SimInfo.MaxV];
    case 'CurrentLayout' % current through junctions
        currs=triu(Sim.Data.Currents{IndexTime});
        Imat=full(abs(currs));
        Cx=Layout.CX(Adj~=0);
        Cy=Layout.CY(Adj~=0);
        Ilist=Imat(Adj~=0);
        I=linspace(0,Sim.SimInfo.MaxI,10*length(Ilist));
        cmap=jet(10*length(Ilist));
        c=interp1(I,cmap,full(Ilist));
        PlotNetworkAux(currAx,X,Y,Cx,Cy,'curr',c);
        clim=[Sim.SimInfo.MinI Sim.SimInfo.MaxI];
    case 'SwitchLayout'  %junction widths
        Wmat=triu(Sim.Data.Wmat{IndexTime});
        Cx=Layout.CX(Adj~=0);
        Cy=Layout.CY(Adj~=0);
        Wlist=Wmat(Adj~=0);
        W=linspace(Sim.SimInfo.MinW,Sim.SimInfo.MaxW,10*length(Wlist));
        cmap=hot(10*length(Wlist));
        c=interp1(W,cmap,full(Wlist));
        PlotNetworkAux(currAx,X,Y,Cx,Cy,'curr',c);
        clim=[Sim.SimInfo.MinW Sim.SimInfo.MaxW];
    case 'ResistanceLayout'
        Rmat=triu(Sim.Data.Rmat{IndexTime});
        Cx=Layout.CX(Adj~=0);
        Cy=Layout.CY(Adj~=0);
        Rlist=Rmat(Adj~=0);
        R=linspace(min([Sim.Settings.Roff Sim.Settings.Ron]),max([Sim.Settings.Roff Sim.Settings.Ron]),10*length(Rlist));
        cmap=flipud(gray(10*length(Rlist)));
        c=interp1(R,cmap,full(Rlist));
        PlotNetworkAux(currAx,X,Y,Cx,Cy,'curr',c);
        clim=[min([Sim.Settings.Roff Sim.Settings.Ron]) max([Sim.Settings.Roff Sim.Settings.Ron])];
        
end
%plot colobar (if necessary) and update its limits
if isempty(findobj(currAx,'Type','ColorBar')) && ~any(strcmp(type,{'ShortestPath','SelEdges'}))
    colormap(currAx,cmap);
    colorbar(currAx);
    caxis(currAx,clim);
end
end



