function PlotGraphLayout(currAx,NodeList,Sim,type,Status,IndexTime,varargin)


glay=varargin{1};
Layout=Sim.SelLayout;
G=Layout.SelGraph;
Adj=(Layout.AdjMat); %adj mat of the simulation

if ~isequal(Status,'redraw')
    
    p=plot(currAx,G,'Layout',glay);
    p.NodeColor='red';
    p.EdgeColor='white';
    p.NodeLabel={};
    switch type
        
        case 'GLayout'
            return;
        case 'GPath'
            return;
        case 'GSTree'
            [T,~] = minspantree(G);
            highlight(p,T,'EdgeColor','r','LineWidth',1.5);
        case 'GDegree'
            G.Nodes.NodeColors = degree(G);
            p.NodeCData = G.Nodes.NodeColors;
            colorbar;
        case 'GVoltage'
            vlist=Sim.Data.Voltages{IndexTime};
            p.NodeCData=full(vlist);
            p.MarkerSize=3;
            colormap(currAx,hot);
            colorbar(currAx);
            caxis([Sim.SimInfo.MinV Sim.SimInfo.MaxV]);
        case 'GCurrent'
            p.MarkerSize=1.5;
            p.LineWidth=1.5;
            currs=(abs(Sim.Data.Currents{IndexTime}));
            [j,i,~]=find(tril(Adj));
            cc=zeros(1,length(j));
            for k=1:length(j)
                cc(k)=currs(i(k),j(k));
            end
            clim=[Sim.SimInfo.MinI Sim.SimInfo.MaxI];
            p.EdgeCData=cc;
            colormap(currAx,gcurrmap);%gcurrmap
            colorbar(currAx);
            caxis(currAx,clim);
        case 'GWidth'
            p.MarkerSize=0.5;
            p.LineWidth=1.5;
            widths=(Sim.Data.Wmat{IndexTime});
            [j,i,~]=find(tril(Adj));
            wd=zeros(1,length(j));
            for k=1:length(j)
                wd(k)=widths(i(k),j(k));
            end
            clim=[Sim.SimInfo.MinW Sim.SimInfo.MaxW];
            p.EdgeCData=wd;
            colormap(currAx,hot);
            colorbar(currAx);
            caxis(currAx,clim);
        case 'GResistance'
            p.MarkerSize=0.5;
            p.LineWidth=1.5;
            res=(Sim.Data.Rmat{IndexTime});
            [j,i,~]=find(tril(Adj));
            wd=zeros(1,length(j));
            for k=1:length(j)
                wd(k)=res(i(k),j(k));
            end
            clim=[min([Sim.Settings.Roff Sim.Settings.Ron]) max([Sim.Settings.Roff Sim.Settings.Ron])];
            p.EdgeCData=wd;
            colormap(currAx,flipud(gcurrmap));%flipud(gray);    
            colorbar(currAx);
            caxis(currAx,clim);
        case 'GBetween'
            c=centrality(G,'betweenness');
            n = numnodes(G);
            p.NodeCData = 2*c./((n-2)*(n-1));
            p.MarkerSize=2;
            colormap(currAx,flip(jet,1));
            colorbar(currAx);
            
    end
else
    p=currAx.Children; % is a graphplt type
    hold(currAx,'on');
    switch type
        case 'GVoltage'
            vlist=Sim.Data.Voltages{IndexTime};
            p.NodeCData=full(vlist);
            p.MarkerSize=3;
        case 'GPath'
            NumEl=height(Sim.Electrodes);
            IdxEl=Sim.Electrodes.PosIndex(NodeList.Value(NodeList.Value<=NumEl));
            if length(IdxEl)>=2
                combs=nchoosek(IdxEl,2); %combs between pair of electrodes
                path=[];
                for i=1:size(combs,1)
                    path=[path shortestpath(G,combs(i,1),combs(i,2))]; %#ok<AGROW>
                end
                %path=unique(path);
                highlight(p,path,'EdgeColor','y','LineWidth',4);
                highlight(p,IdxEl,'NodeColor','g','MarkerSize',3)
                
            end
        case 'GMaxflow'
            p.EdgeColor='white';
            p.LineWidth=0.5;
            NumEl=height(Sim.Electrodes);
            IdxEl=Sim.Electrodes.PosIndex(NodeList.Value(NodeList.Value<=NumEl));
            if length(IdxEl)==2
                [~,path]=maxflow(G,IdxEl(1),IdxEl(2));
                highlight(p,path,'EdgeColor','g','LineWidth',4);
                highlight(p,IdxEl,'NodeColor','g','MarkerSize',3)
            end
            
            
        case 'GCurrent'
            p.EdgeColor='none';
            currs=(abs(Sim.Data.Currents{IndexTime}));
            [j,i,~]=find(tril(Adj));
            cc=zeros(1,length(j));
            for k=1:length(j)
                cc(k)=currs(i(k),j(k));
            end
            p.EdgeCData=cc;
        case 'GWidth'
            p.EdgeColor='none';
            p.MarkerSize=0.5;
            p.LineWidth=1.5;
            widths=(Sim.Data.Wmat{IndexTime});
            [j,i,~]=find(tril(Adj));
            wd=zeros(1,length(j));
            for k=1:length(j)
                wd(k)=widths(i(k),j(k));
            end
            p.EdgeCData=wd;
        case 'GResistance'
            p.EdgeColor='none';
            p.MarkerSize=0.5;
            p.LineWidth=1.5;
            res=(Sim.Data.Rmat{IndexTime});
            [j,i,~]=find(tril(Adj));
            wd=zeros(1,length(j));
            for k=1:length(j)
                wd(k)=res(i(k),j(k));
            end
            p.EdgeCData=wd;
            
            
            
            
    end
    
end

%% plot electrode markers
switch type
    case {'GWidth','GResistance','GLayout','GCurrent'}
        DSel=contains(Sim.Electrodes.Name,'Drain');
        SSel=contains(Sim.Electrodes.Name,'Source');
        highlight(p,Sim.Electrodes.PosIndex(DSel),'NodeColor','y','MarkerSize',6);
        labelnode(p,Sim.Electrodes.PosIndex(DSel),Sim.Electrodes.Name(DSel));
        highlight(p,Sim.Electrodes.PosIndex(SSel),'NodeColor','g','MarkerSize',6);
        labelnode(p,Sim.Electrodes.PosIndex(SSel),Sim.Electrodes.Name(SSel));
end



end