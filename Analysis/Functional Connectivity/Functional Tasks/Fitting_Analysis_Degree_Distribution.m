cd('C:\Users\61424\Documents\GitHub\CODE\Data\Functional Connectivity')
load('Degree.mat')

for i=1:length(values)
    r = double(values{i});
    
    gamfitVals{i}=gamfit(r);
    expfitVals{i}=expfit(r);
    wblfitVals{i}=wblfit(r);
    lognfitVals{i}=lognfit(r);%phat
    
    [logLg,ncoV{i}] = gamlike(gamfit(r),r);
    logLe{i} = explike(expfit(r),r);
    logLw = wbllike(wblfit(r),r);
    [logLn,pcoV{i}] = lognlike(lognfit(r),r);     
    [~,~,logLp] = plfit(r);
    
    [~,bicg(i)] = aicbic(logLg,1:1:size(logLg,1),length(r));
    [~,bice(i)] = aicbic(logLe{i},1:1:size(logLe{i},1),length(r));
    [~,bicw(i)] = aicbic(logLw,1:1:size(logLw,1),length(r));
    [~,bicp(i)] = aicbic(-1*logLp,1:1:size(logLp,1),length(r));
    [~,bicn(i)] = aicbic(logLn,1:1:size(logLn,1),length(r));

    
    if bicg(i) < bice(i) & bicg(i) < bicw(i) & bicg(i) < bicp(i) & bicg(i) < bicn(i)
        winner{i}='Gamma';%sprintf('%s%d','Gamma wins! BIC = ',bicg)
    elseif bice(i) < bicg(i) & bice(i) < bicw(i) & bice(i) < bicp(i) & bice(i) < bicn(i)
        winner{i}='Exponential';%sprintf('%s%d','Exponential wins! BIC = ',bicg)
    elseif bicw(i) < bicg(i) & bicw(i) < bice(i) & bicw(i) < bicp(i) & bicw(i) < bicn(i)
        winner{i}='Stretched Exponential';%sprintf('%s%d','Stretched Exponential wins! BIC = ',bicg)
    elseif bicp(i) < bicg(i) & bicp(i) < bicw(i) & bicp(i) < bicp(i) & bicp(i) < bicn(i)
        winner{i}='Power Law';%sprintf('%s%d','Power Law wins! BIC = ',bicg)
     elseif bicn(i) < bicg(i) & bicn(i) < bicw(i) & bicn(i) < bicp(i) & bicn(i) < bice(i)
        winner{i}='Log Normal';%sprintf('%s%d','Power Law wins! BIC = ',bicg)
    end
end 
%%  
% Overlay plot fit for everything:
for i =1%:length(values)
%    
muhat=lognfitVals{i}(1);
sigmahat=lognfitVals{i}(2);
figure()
%    h{i}=histogram(values{i})
   histfit(values{i},20,'normal')
   hold on
   h2{i}=histfit(double(values{i}),20,'gamma')
   delete(h2{i}(1))
   h2{i}(2).Color ='g';
   h2{i}=histfit(double(values{i}),20,'exponential')
   delete(h2{i}(1))
   h2{i}(2).Color ='m';
   h2{i}=histfit(double(values{i}),20,'logn')
   delete(h2{i}(1))
   h2{i}(2).Color ='k';
   legend('Degree','Normal Fit','Gamma Fit','Exponential','Log Normal')
%    hold on
%    xx = linspace(0,double(max(values{i})));
%    line(xx,gamcdf(xx,gamfitVals{i}(1),gamfitVals{i}(2),ncoV{i}),'color','r')
%    line(xx,expcdf(xx,expfitVals{i}(1),logLe{i}),'color','g')
%    line(xx,logncdf(xx,muhat,sigmahat,pcoV{i}),'color','m')
end 
