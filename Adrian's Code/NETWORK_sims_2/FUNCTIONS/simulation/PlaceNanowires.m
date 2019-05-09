function NewNet=PlaceNanowires(InfoText,Net)
%% generate nanowires and place them randomly in space

rng('shuffle');

X=Net.SizeX*rand(Net.Number,1);
Y=Net.SizeY*rand(Net.Number,1);
if isequal(Net.AngleCheck,0)
    Angles=-Net.Angles+(2*Net.Angles)*rand(Net.Number,1);
else
    Angles=-Net.Angles+(2*Net.Angles)*round(rand(Net.Number,1));
end
theta=deg2rad(Angles);
Lengths=Net.Length*Net.Disp/100*randn(1,length(X))...
            +Net.Length;
x1=zeros(length(X),1);
x2=x1;
y1=x1;
y2=x1;
if isequal(Net.ManCheck,0)
    for i=1:length(X)
        x1(i)=X(i);
        x2(i)=X(i)+Lengths(i)*cos(theta(i));
        y1(i)=Y(i);
        y2(i)=Y(i)+Lengths(i)*sin(theta(i));
    end
else
    uiwait(msgbox(strcat('Manually place beggining and ending position of  ',num2str(length(X)), ' nanowires')));
    set(gca,'XLim',[0 Net.SizeX]);
    set(gca,'YLim',[0 Net.SizeY]);
    hold all
    for i=1:length(X)        
        [xb,yb]=ginput(1);
        scatter(xb,yb,20,'r');
        [xe,ye]=ginput(1);
        scatter(xe,ye,20,'r');
        plot([xb,xe],[yb,ye]);
        InfoText.String=strcat(num2str(length(X)-i),' Nanowires yet to place, hurry up!');
        x1(i)=xb;x2(i)=xe;
        y1(i)=yb;y2(i)=ye;
        
    end
    
end
InfoText.String='Network completed';

%% ordering nanowires by initial x position
[x1,k]=sort(x1);
x2=x2(k);
y1=y1(k);
y2=y2(k);

%%saving positionss as table
NewNet=table(x1,x2,y1,y2);

end