function SimInfo=GetSimInfo(SimData)

%this could be done easily with cellfun,
% but sparrays are not supported by it

Lc=length(SimData.Currents);
minc=zeros(Lc,1);
maxc=zeros(Lc,1);
minw=zeros(Lc,1);
maxw=zeros(Lc,1);
minv=zeros(Lc,1);
maxv=zeros(Lc,1);
for i=1:Lc
    [~,~,c]=find(SimData.Currents{i});
    [~,~,w]=find(SimData.Wmat{i});
    [~,~,v]=find(SimData.Voltages{i});
    
    if isempty(c)
        minc(i)=0;
        maxc(i)=0;
    else
        minc(i)=min(abs(c));
        maxc(i)=max(abs(c));
    end
    minw(i)=min(w);
    maxw(i)=max(w);
    if isempty(v)
        minv(i)=0;
        maxv(i)=0;
    else
        minv(i)=min(v);
        
        maxv(i)=max(v);
    end
    
end
SimInfo.MaxI=max(maxc);
SimInfo.MinI=min(minc);
SimInfo.MaxW=max(maxw);
SimInfo.MinW=min(minw);
SimInfo.MinV=min(minv);
SimInfo.MaxV=max(maxv);

end