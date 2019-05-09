function PlotNetwork(NetStruct,choices)

cla;
Network=NetStruct.LayOut;

x1=Network.x1;
x2=Network.x2;
y1=Network.y1;
y2=Network.y2;
Cx=Network.XInt; Cx=Cx(~cellfun(@isempty,Cx));
Cy=Network.YInt; Cy=Cy(~cellfun(@isempty,Cy));
Cx=cell2num(Cx);
Cy=cell2num(Cy);
currAx=gca;
X=[x1';x2'];
Y=[y1';y2'];
NoW=length(x1);




if isequal(choices,'domains')
    Domains=NetStruct.Domains;
    LDom=length(Domains);
    
    %%color map choices should be here
    cmap=lines(LDom);
    cm=zeros(NoW,3);
    plot(currAx,X,Y,'b','LineWidth',2);
    for i=1:LDom
        idx=Domains{i};
        cm(idx,:)=repmat(cmap(i,:),[length(idx) 1]);
    end
    cm=flipud(cm);
    set(currAx.Children,{'Color'},num2cell(cm,2));
else
    plot(currAx,X,Y,'b');
end
hold on
scatter(Cx,Cy,2,'r');



end