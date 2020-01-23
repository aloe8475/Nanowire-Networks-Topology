%cll
global glo;
%  rqa_stat - RQA statistics - [recrate DET LMAX ENT TND LAM TT]
contactn=50;
glo.a=0;glo.plotpsd=0;glo.dutyratio=0;glo.ap=0;glo.c=0;
glo.b=0;glo.M=0;glo.t=0;glo.testcurrent=0;glo.contactn=contactn;
f=67;
N=eye(10);N=[N;N];
% for dutyratio=20:10:70
%     a=4./dutyratio;
% MonolithicDemoKevin;
% % Output.networkResistance(1:500)=[];
% % rps(Output.networkResistance);
% end
% % SimulationOptions.T.*1./2
glo.a=1;
glo.dutyratio=50;
glo.b=0.015;
glo.th=0.5;
glo.c=70;
glo.f=149;
glo.ap=1;
glo.tau=5;
glo.testn=40;
% i=1;
%  MonolithicDemoKevin;
% for a=2
% glo.contactn=a;
% MonolithicDemoKevin;
% mse(i)=Output.mse;
% i=i+1;
% end
% mse=zeros(1,6);
% for  tau=7
%
% %     tau=10;
% MonolithicDemoKevin(glo);
% mse(i)=glo.contactn;
% i=i+1;
% end
% end
if 0
    MonolithicDemoKevin;
    recurent=Output.networkCurrent(end,:);recurent(:,2)=[];
    cc=ann54(recurent);
    [maxVal maxInd] =max(cc);
    bb(b)=maxInd;
end
if 0
    for b=1:10
        MonolithicDemoKevin;
        recurent=Output.networkCurrent(end,:);
        re=recurent;re(:,2)=[];
        cc=ann54(re);
        [maxVal maxInd] =max(cc);
        bb(b)=maxInd;
    end
end
if 0
    for b=1:10
        MonolithicDemoKevin;
        recurent(b,:)=Output.networkCurrent(end,:);
    end
    recurent(:,2)=[];
    recurent=[recurent;recurent];
    % save dat1 recurent;
    % save dat2 N;
end
% figure
% x=Connectivity.VertexPosition(:,1);y=Connectivity.VertexPosition(:,2);
% n=length(x);
% d = linspace(1,10,length(x));
% scatter(x,y,'filled','o');
% text(x,y,arrayfun(@(x)['  ' num2str(x)],1:n,'UniformOutput',0));
%
% MonolithicDemoKevinMultiElec
%    MonolithicDemoKevin;
if 0
    for ap=0.7:0.05:0.7
        glo.ap=ap;
        MonolithicDemoKevin;
        % record every junctions' V and G
        absV=Output.wireVoltage';
        for i=1:length(snapshots)
            absV(:,i)=getAbsoluteVoltage(snapshots{floor(i)},Connectivity,SimulationOptions.ContactNodes);
        end
        recordjvoltage=zeros(length(snapshots{1}.Voltage),length(snapshots));
        for i=1:length(snapshots)
            recordjvoltage(:,i)=snapshots{i}.Voltage;
        end
        recordjG=zeros(length(snapshots{1}.Resistance),length(snapshots));
        for i=1:length(snapshots)
            recordjG(:,i)=1./snapshots{i}.Resistance;
        end
        afplot(f,recordjvoltage,recordjG,snapshots{floor(1)},Stimulus,absV);
        figure
        x=Connectivity.EdgePosition(:,1);y=Connectivity.EdgePosition(:,2);
        n=length(x);
        d = linspace(1,10,length(x));
        scatter(x,y,'filled','o');
        text(x,y,arrayfun(@(x)['  ' num2str(x)],1:n,'UniformOutput',0));
        figure
        x=Connectivity.VertexPosition(:,1);y=Connectivity.VertexPosition(:,2);
        n=length(x);
        d = linspace(1,10,length(x));
        scatter(x,y,'filled','o');
        text(x,y,arrayfun(@(x)['  ' num2str(x)],1:n,'UniformOutput',0));
        % close all;
    end
end
if 1
    MonolithicDemoKevin;
    absV=Output.wireVoltage';
    %     MonolithicDemoKevin;
    % linear regression part
    % target=Stimulus.Signal;help
    stepm=[1000,500,250,200,125,100,50,10];
    target= glo.ap*square(2*glo.a*pi*Stimulus.Frequency*(Stimulus.TimeAxis));
    [mse,accuracy]=multipleregression(target,absV,stepm,Stimulus);
    xx=1000./stepm(1:8);
    figure, plot(xx,accuracy,'o-');
    xlabel('times of regression');
    ylabel('accuracy','rotation',0);
    set(gca,'FontSize',16);  set(gca,'XMinorTick','on','YMinorTick','on');
end
%
% %  system('shutdown -s');
%
% % rps(Output.networkResistance);
%
% figure,plot(Stimulus.TimeAxis,Output.networkCurrent.*Stimulus.Signal);
% figure;
% plot(Output.networkCurrent,Stimulus.Signal);?
% end
% dutyratio=dutyratio+4;
% end
% Output.networkResistance(1:300)=[];
%  rps(Output.networkResistance);