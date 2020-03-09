function sclist=cross_correlator(avelist,targetlist)
sclist=zeros(1,length(avelist));
for i=1:length(avelist)
    cave=avelist{i};
    ctar=NumbToBoolArray(targetlist(i));
    
    un=(cave-mean(cave));
    dos=(ctar-mean(ctar));
    
    u=sum(un.*dos);
    d=sqrt(sum(un.^2)*sum(dos.^2));
    
    sclist(i)=u/d;
end