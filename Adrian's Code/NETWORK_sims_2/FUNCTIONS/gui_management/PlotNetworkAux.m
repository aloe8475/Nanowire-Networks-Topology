function PlotNetworkAux(currAx,X,Y,Cx,Cy,type,varargin)
switch type
    case 'all'
        
        plot(currAx,X,Y,'b');
        hold on
        scatter(currAx,Cx,Cy,2,'r');
    case 'electrodes'
        Annotate=varargin{2};
        
        lc=findobj(currAx,'LineWidth',2.2);
        tx=findobj(currAx,'Type','text');
        if ~isempty(lc)
            delete(lc);
            
        end
        if ~isempty(tx)
            delete(findobj(tx,'Color','green'));
        end
        plot(currAx,X,Y,'yellow','LineWidth',2.2);
        if ~isequal(Annotate,0)
            ElectrodeText=varargin{1};
            text(currAx,Cx+0.7,Cy+0.7,ElectrodeText,'Color','green','FontSize',8,...
                'Interpreter','none','FontWeight','bold');
        end
        
    case 'nodes'
        
        Annotate=varargin{2};
        
        sc=findobj(currAx,'MarkerFaceColor','white');
        tx=findobj(currAx,'Type','text');
        if ~isempty(sc)
            delete(sc);
        end
        if ~isempty(tx)
            delete(findobj(tx,'Color','white'));
        end
        scatter(currAx,Cx,Cy,50,'w','h','MarkerFaceColor','white');
        if ~isequal(Annotate,0)
            NodeText=varargin{1};
            NodeText=strrep(NodeText,'Node_','n');
            text(currAx,Cx+0.4,Cy-0.3,NodeText,'Color','white','FontSize',7,...
                'Interpreter','none','FontWeight','bold');
            
        end
        
    case 'lines'
        plot(currAx,X,Y,'b');
        
    case 'crosses'
        scatter(currAx,Cx,Cy,2,'r');
        
    case 'path'
        color=varargin{1};
        plot(currAx,X,Y,color,'LineWidth',3);
    case 'volt'
        color=varargin{1};
        p=plot(currAx,X,Y,'b','LineWidth',3);
        %color=flipud(color);
        set(p,{'Color'},num2cell(color,2));
    case 'curr'
        color=varargin{1};
        sc=scatter(currAx,Cx,Cy,25,'s');
        sc.MarkerEdgeColor='flat';
        sc.MarkerFaceColor='flat';
        sc.CData=color;
end