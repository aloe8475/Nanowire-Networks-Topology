function varargout=fillRfromIV(condition,varargin)
%% fill R from calculated IV after kirchoff solver.
% this takes any kind of first order time evolution
% to generate evolution of R following memristor equation



switch condition
    case 'ini'
        sims=varargin{1};
        settings=sims.Settings;        
        simout=table();
        simout.AdjMat{1}=sims.SelLayout.AdjMat;
        [i,j,~]=find(simout.AdjMat{1});
        if isempty(sims.LastW)
        ini_widths=settings.IniW*settings.SigmaW*randn(1,length(i))...
            +settings.IniW;
        newW=sparse(i,j,ini_widths);
        newW(newW>=settings.MaxW)=settings.MaxW;
        newW(newW<=0)=0;   
        newW=triu(newW)'+triu(newW);
        else
            newW=sims.LastW;
        end
        simout.Rmat{1}=sparse((settings.Ron-settings.Roff)/settings.MaxW.*newW+settings.Roff.*simout.AdjMat{1});
        
        [i,j,r]=find(simout.Rmat{1});
        g=1./r;
        simout.Gmat{1}=sparse(i,j,g);
        
        simout.Wmat{1}=newW;
        
       
        varargout{1}=simout;
        
        
    case 'loop'
        sim=varargin{1};settings=varargin{2};
        newsim=UpdateRmat(sim,settings);
        varargout{1}=newsim;

    otherwise
        return;
end



end