function sclist=Euclid_distance(avelist,targetlist)
sclist=zeros(1,length(avelist));
for i=1:length(avelist)
    cave=avelist{i};
    cave=minMax(cave,0.01);
    ctar=NumbToBoolArray(targetlist(i));
    u=fuzzy_bool(cave,'labview');
    d=fuzzy_bool(ctar,'labview');
    sclist(i)=abs(u-d);
end
end