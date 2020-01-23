function [mse,accuracy]=multipleregression(target,absV,stepm,Stimulus)
inputx=[ones(length(target),1),absV']; 
for o=1:6  % o9 marks the times of reression depends on stepm
i=1;step=stepm(o); 
for ii=1:step:1/(Stimulus.Frequency*Stimulus.dt)%for i in range(times_of_regresssion), ii marks the begin point of separated inputx, thus ii = (i-1)*step+1 
    recordv=[]; targetii=[];%redo this 2 mid-matrix everytime
        for iii=0:floor(length(target)/(1/(Stimulus.Frequency*Stimulus.dt)))-1
            targetii=[targetii; target(ii+iii*1/(Stimulus.Frequency*Stimulus.dt):ii+iii*1/(Stimulus.Frequency*Stimulus.dt)+step-1)];
            recordv=[recordv,absV(:,ii+iii*1/(Stimulus.Frequency*Stimulus.dt):ii+iii*1/(Stimulus.Frequency*Stimulus.dt)+step-1)];
        end%I am collecting every period of signal to do regression, not just one, like there are 10 period, iii=0:9
        inputxii=[ones(length(targetii),1),recordv'];
        [a1(:,i),bint,r,rint,state]=regress(targetii,inputxii);%regression, noticed that a1 have i column
        i=i+1;
end
a1=[a1(:,end) a1(:,1:end-1)];
% load dat a1
for i=1:length(Stimulus.TimeAxis)
result(i)=inputx(i,:)*a1(:,rem(ceil(i/step),1/(Stimulus.Frequency*Stimulus.dt*step))+1); % every part(i part) doing regression independently
end
% figure
% subplot(2,1,1),plot(target); subplot(2,1,2),plot(result);set(gca,'FontSize',16);  set(gca,'XMinorTick','on','YMinorTick','on');
figure,%plot the result
plot(target,'r');
hold on,
plot(result,'b--');ylim([-1.5,1.5]);set(gca,'FontSize',16);  set(gca,'XMinorTick','on','YMinorTick','on');
xlabel('t');
ylabel('v','rotation',0);
hold off
mse(o) = sum((target-result').^2)./length(target);
accuracy(o)=1-sqrt(sum((target-result').^2)./sum(target.^2));
% nrmse(o) = sqrt(sum((1-result.\target').^2)/length(target));
% accuracy(o) =1-mse(o);
% disp( ['MSE = ', num2str( mse )] );disp( ['accuracy = ', num2str( accuracy )] );
end
end