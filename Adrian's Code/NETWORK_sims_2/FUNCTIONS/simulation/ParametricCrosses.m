function NewNet=ParametricCrosses(Network)

%% Resolve system of equations of one nanowire with his suspected
%% neighbors, if the parametric solution falls within [0 1], then the 
%% two nanowires intersect

x1=Network.x1;
x2=Network.x2;
y1=Network.y1;
y2=Network.y2;
CrossIndex=Network.CrossIndex;
NoW=length(x1);
Network.XInt=cell(NoW,1);
Network.YInt=cell(NoW,1);
Network.IndexInt=cell(NoW,1);
Network.NEdges=zeros(NoW,1);
for i=1:NoW
%%select candidates around nanowire i    
Set=CrossIndex{i};

% get origin, length
p0=[x1(i);y1(i)];
s0=[x2(i)-x1(i);y2(i)-y1(i)]; 

p=[x1(Set)'; y1(Set)'];
s=[(x2(Set)-x1(Set))'; (y2(Set)-y1(Set))'];

%initialize matrices    
ls=length(s(1,:));
A=zeros(ls*2,1);
B=zeros(ls*2,ls*2);
p0x=p0(1);p0y=p0(2);
s0x=s0(1);s0y=s0(2);
px=p(1,:);py=p(2,:);
sx=s(1,:);sy=s(2,:);

%fill parametric matrices
for j=1:(ls)
    A(2*j-1)=p0x-px(j);
    A(2*j)=p0y-py(j);
   
    B(2*j-1,2*j-1)=-s0x;B(2*j-1,2*j)=sx(j);
    B(2*j,2*j-1)=-s0y;B(2*j,2*j)=sy(j);
end
%solve system of equations
sols=B\A;

%fill succesful candidates and store in network table
n=0;
XIntersect=[];YIntersect=[];IndexIntersect=[];
for j=1:ls
    t=sols(2*j-1);
    if t>=0 && t<=1
        s=sols(2*j);
        if s>=0 && s<=1
           n=n+1;
           XIntersect(n)=p0x+t*s0x; %#ok<*AGROW>
           YIntersect(n)=p0y+t*s0y;
           IndexIntersect(n)=Set(j);
         
        end
    end
end

if ~isempty(IndexIntersect)
    Network.XInt(i)={XIntersect};
    Network.YInt(i)={YIntersect};
    Network.IndexInt(i)={IndexIntersect};
    Network.NEdges(i)=length(IndexIntersect);
end
end
Network.CrossIndex=[];
NewNet=Network;
end