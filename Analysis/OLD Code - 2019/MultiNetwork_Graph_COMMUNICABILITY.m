load('C:\Users\aloe8475\Documents\PhD\GitHub\CODE\Data\Explore Analysis\Ordered+Random Data\WScomm_curr_plots.mat');

dbstop if error


clear randNoise

save_directory='C:\Users\aloe8475\Documents\PhD\GitHub\CODE\Data\Figures\Explore Analysis\Cross-Network Explore\Graph Theory\';
%COLORS
AgNW.s100.length=length(Curr{2});
AgNW.s100.color_1=hex2rgb('00ffff');
AgNW.s100.color_end=hex2rgb('0049E9');
AgNW.s100.colors=[linspace(AgNW.s100.color_1(1),AgNW.s100.color_end(1),AgNW.s100.length)', linspace(AgNW.s100.color_1(2),AgNW.s100.color_end(2),AgNW.s100.length)', linspace(AgNW.s100.color_1(3),AgNW.s100.color_end(3),AgNW.s100.length)'];

AgNW.s300.length=length(Curr{3});
AgNW.s300.color_1=hex2rgb('ffff84');
AgNW.s300.color_end=hex2rgb('FC9C03');
AgNW.s300.colors=[linspace(AgNW.s300.color_1(1),AgNW.s300.color_end(1),AgNW.s300.length)', linspace(AgNW.s300.color_1(2),AgNW.s300.color_end(2),AgNW.s300.length)', linspace(AgNW.s300.color_1(3),AgNW.s300.color_end(3),AgNW.s300.length)'];

AgNW.s500.length=length(Curr{4});
AgNW.s500.color_1=hex2rgb('01ff07');
AgNW.s500.color_end=hex2rgb('00551A');
AgNW.s500.colors=[linspace(AgNW.s500.color_1(1),AgNW.s500.color_end(1),AgNW.s500.length)', linspace(AgNW.s500.color_1(2),AgNW.s500.color_end(2),AgNW.s500.length)', linspace(AgNW.s500.color_1(3),AgNW.s500.color_end(3),AgNW.s500.length)'];


WS.b1.length=length(Curr{1})/3;
WS.b05.length=length(Curr{1})/3;
WS.b015.length=length(Curr{1})/3;

WS.b1.color_1=hex2rgb('FFE8B6');
WS.b1.color_last=hex2rgb('FFAE00');
WS.b1.colors = [linspace(WS.b1.color_1(1),WS.b1.color_last(1),WS.b1.length)', linspace(WS.b1.color_1(2),WS.b1.color_last(2),WS.b1.length)', linspace(WS.b1.color_1(3),WS.b1.color_last(3),WS.b1.length)'];

WS.b05.color_1=hex2rgb('FFDFBD');
WS.b05.color_last=hex2rgb('FF8300');
WS.b05.colors = [linspace(WS.b05.color_1(1),WS.b05.color_last(1),WS.b05.length)', linspace(WS.b05.color_1(2),WS.b05.color_last(2),WS.b05.length)', linspace(WS.b05.color_1(3),WS.b05.color_last(3),WS.b05.length)'];

WS.b015.color_1=hex2rgb('FFCCB0');
WS.b015.color_last=hex2rgb('FF5900');
WS.b015.colors = [linspace(WS.b015.color_1(1),WS.b015.color_last(1),WS.b015.length)', linspace(WS.b015.color_1(2),WS.b015.color_last(2),WS.b015.length)', linspace(WS.b015.color_1(3),WS.b015.color_last(3),WS.b015.length)'];

cElegans.color=hex2rgb('FC9C03');

%DATA
AgNW.s100.data.current=Curr{2};
AgNW.s300.data.current=Curr{3};
AgNW.s500.data.current=Curr{4};
cElegans.data.current=Curr{5};
WS.b015.data.current=Curr{1}(1:9);
WS.b05.data.current=Curr{1}(10:18);
WS.b1.data.current=Curr{1}(19:27);
WS.b015.data.comm=COMM{1}(1:9);
WS.b05.data.comm=COMM{1}(10:18);
WS.b1.data.comm=COMM{1}(19:27);
AgNW.s100.data.comm=COMM{2};
AgNW.s300.data.comm=COMM{3};
AgNW.s500.data.comm=COMM{4};
cElegans.data.comm=COMM{5};

for i = 1:7
    switch i
        case 1
            for j = 1:length(WS.b015.data.comm)
                WS.b015.data.comm{j}=WS.b015.data.comm{j}(WS.b015.data.comm{j}~=0);
                randNoise{i,j}=0.9+0.2*rand(1,length(WS.b015.data.comm{j})); %create random noise
            end
        case 2
            for j = 1:length(WS.b05.data.comm)
                WS.b05.data.comm{j}=WS.b05.data.comm{j}(WS.b05.data.comm{j}~=0);
                randNoise{i,j}=1.9+0.2*rand(1,length(WS.b05.data.comm{j})); %create random noise
            end
        case 3
            for j = 1:length(WS.b1.data.comm)
                WS.b1.data.comm{j}=WS.b1.data.comm{j}(WS.b1.data.comm{j}~=0);
                
                randNoise{i,j}=2.9+0.2*rand(1,length(WS.b1.data.comm{j})); %create random noise
            end
        case 4
            for j = 1:length(AgNW.s100.data.comm)
                AgNW.s100.data.comm{j}=AgNW.s100.data.comm{j}(AgNW.s100.data.comm{j}~=0);
                randNoise{i,j}=3.9+0.2*rand(1,length(AgNW.s100.data.comm{j})); %create random noise
            end
        case 5
            for j = 1:length(AgNW.s300.data.comm)
                AgNW.s300.data.comm{j}=AgNW.s300.data.comm{j}(AgNW.s300.data.comm{j}~=0);
                
                randNoise{i,j}=4.9+0.2*rand(1,length(AgNW.s300.data.comm{j})); %create random noise
            end
            
        case 6
            for j = 1:length(AgNW.s500.data.comm)
                AgNW.s500.data.comm{j}=AgNW.s500.data.comm{j}(AgNW.s500.data.comm{j}~=0);
                
                randNoise{i,j}=5.9+0.2*rand(1,length(AgNW.s500.data.comm{j})); %create random noise
            end
            
        case 7
            cElegans.data.comm=cElegans.data.comm(cElegans.data.comm~=0);
            randNoiseCEleg=6.9+0.2*rand(1,length(cElegans.data.comm)); %create random noise
    end
end


%% Mean comparisons:

% %All network densities:
%  FigList = allchild(groot);
%         categories.COMM={vertcat(AgNW.s100.data.comm{:})  vertcat(AgNW.s300.data.comm{:}) vertcat(WS.b015.data.comm{:}) vertcat(WS.b05.data.comm{:}) vertcat(WS.b1.data.comm{:}) cElegans.data.comm};
%         maxlength = max(cellfun(@numel, categories.COMM));
%         COMMtemp = cellfun(@(v) [v', nan(1, maxlength-numel(v))], categories.COMM, 'UniformOutput', false);
% 
% [ANOVA.COMM.p,ANOVA.COMM.AnovaTab,ANOVA.COMM.Stats]=anova1([COMMtemp{1}' COMMtemp{2}' COMMtemp{3}' COMMtemp{4}' COMMtemp{5}' COMMtemp{6}'],{'100nw AgNW','300nw AgNW','WS {\beta} = 0.15','WS {\beta} = 0.5','WS {\beta} = 1', '{\it C. Elegans}'});
%     [POSTHOC.COMM.c, POSTHOC.COMM.m, POSTHOC.COMM.h, POSTHOC.COMM.nms]=multcompare(ANOVA.COMM.Stats,'alpha',.05/6,'ctype','bonferroni');
  

%%
%% COMM VS CURR:
f = figure('Position',[0 0 1920 1080]);

% for i = 1:length(Curr)
    for j = 1:length(Curr{2})
%         if i == 1
%             if j < 10
%                 g1{j}=gscatter([WS.b015.data.current{j}],([WS.b015.data.comm{j}]),[],AgNW.s300.colors(j,:),[],15,'off');
%                     set(gca, 'YScale', 'log')
% 
%                 hold on;
%             elseif j < 19
%                 g2{j}=gscatter([WS.b05.data.current{j-9}]',([WS.b05.data.comm{j-9}]),[],AgNW.s300.colors(j-9,:),[],15,'off');
%                     set(gca, 'YScale', 'log')
% 
%                 hold on;
%             else
%                 g3{j}=gscatter([WS.b1.data.current{j-18}]',([WS.b1.data.comm{j-18}]),[],AgNW.s300.colors(j-18,:),[],15,'off');
%                     set(gca, 'YScale', 'log')
% 
%                 hold on;
%             end
%         elseif i == 2
%             current{j}=AgNW.s100.data.current{j}(AgNW.s100.data.current{j}>0);
%             comm{j}=AgNW.s100.data.comm{j}(AgNW.s100.data.current{j}>0);
            g4{j}=gscatter(Curr{2}{j},COMM{2}{j})%,AgNW.s100.colors(j,:),[],15,'off');
            set(g4{j},'MarkerSize',15);
%                 set(gca, 'YScale', 'log')
                set(gca,'XScale','log')
            hold on;
%         elseif i == 3
%             g5{j}=gscatter([AgNW.s300.data.current{j}]',([AgNW.s300.data.comm{j}]),[],AgNW.s300.colors(j,:),[],15,'off');
%                 set(gca, 'YScale', 'log')

%             hold on;
%         elseif i == 4
%             g6{j}=gscatter([AgNW.s500.data.current{j}]',([AgNW.s500.data.comm{j}]),[],AgNW.s500.colors(j,:),[],15,'off');
%                 set(gca, 'YScale', 'log')
% 
%             hold on;
            
%         elseif i == 5
%             g7=gscatter(randNoiseCEleg', (cElegans.data.comm), [],cElegans.color,[],15,'off');
%                 set(gca, 'YScale', 'log')
% 
%             hold on;
%         end
%     end
end
legend;
xlabel('log(COMMCollective)');
ylabel('log(Current)');
% xlim([0.5 7.5])
% yticks([-10:2.5:15])
% xticklabels({"300 Node WS {\beta} = 0.15","300 Node WS {\beta} = 0.5","300 Node WS {\beta} = 1","100nw","300nw","500nw","{\it C. Elegans}"});
% ylabel("log(Communicability)");

% print(f,'-painters','-dpng','-r300',[save_directory 'Communicability WS vs AgNW vs C Elegans.png']);
% print(f,'-painters','-dpdf','-bestfit','-r150',[save_diclrectory 'Communicability WS vs AgNW vs C Elegans.pdf']);

%% PLOTs:
f = figure('Position',[0 0 1920 1080]);

for i = 1:length(Curr)
    for j = 1:length(Curr{i})
        if i == 1
            if j < 10
                g1{j}=gscatter([randNoise{1,j}]',([WS.b015.data.comm{j}]),[],AgNW.s300.colors(j,:),[],15,'off');
                    set(gca, 'YScale', 'log')

                hold on;
            elseif j < 19
                g2{j}=gscatter([randNoise{2,j-9}]',([WS.b05.data.comm{j-9}]),[],AgNW.s300.colors(j-9,:),[],15,'off');
                    set(gca, 'YScale', 'log')

                hold on;
            else
                g3{j}=gscatter([randNoise{3,j-18}]',([WS.b1.data.comm{j-18}]),[],AgNW.s300.colors(j-18,:),[],15,'off');
                    set(gca, 'YScale', 'log')

                hold on;
            end
        elseif i == 2
            g4{j}=gscatter([randNoise{4,j}]',([AgNW.s100.data.comm{j}]),[],AgNW.s100.colors(j,:),[],15,'off');
                set(gca, 'YScale', 'log')

            hold on;
        elseif i == 3
            g5{j}=gscatter([randNoise{5,j}]',([AgNW.s300.data.comm{j}]),[],AgNW.s300.colors(j,:),[],15,'off');
                set(gca, 'YScale', 'log')

            hold on;
        elseif i == 4
            g6{j}=gscatter([randNoise{6,j}]',([AgNW.s500.data.comm{j}]),[],AgNW.s500.colors(j,:),[],15,'off');
                set(gca, 'YScale', 'log')

            hold on;
            
        elseif i == 5
            g7=gscatter(randNoiseCEleg', (cElegans.data.comm), [],cElegans.color,[],15,'off');
                set(gca, 'YScale', 'log')

            hold on;
        end
    end
end


xlim([0.5 7.5])
% yticks([-10:2.5:15])
xticklabels({"300 Node WS {\beta} = 0.15","300 Node WS {\beta} = 0.5","300 Node WS {\beta} = 1","100nw","300nw","500nw","{\it C. Elegans}"});
ylabel("log(Communicability)");

print(f,'-painters','-dpng','-r300',[save_directory 'Communicability WS vs AgNW vs C Elegans.png']);
% print(f,'-painters','-dpdf','-bestfit','-r150',[save_diclrectory 'Communicability WS vs AgNW vs C Elegans.pdf']);


%% SAVE TABLES
% %% SPLIT INTO THE 9-11 DIFFERENT NETWORK SPARSITIES:
% for i = 1:9
%     FigList = allchild(groot);
% %     categories{i}.COMM={vertcat(AgNW.s300.data.comm{i}) vertcat(WS.b015.data.comm{i}) vertcat(WS.b05.data.comm{i}) vertcat(WS.b1.data.comm{i}) cElegans.data.comm};
% %     maxlength = max(cellfun(@numel, categories{i}.COMM));
% %     COMMtemp = cellfun(@(v) [v', nan(1, maxlength-numel(v))], categories{i}.COMM, 'UniformOutput', false);
% %     [ANOVA.COMM{i}.p,ANOVA.COMM{i}.AnovaTab,ANOVA.COMM{i}.Stats]=anova1([COMMtemp{1}' COMMtemp{2}' COMMtemp{3}' COMMtemp{4}' COMMtemp{5}'],{'300nw AgNW','WS beta = 0.15','WS beta = 0.5','WS beta = 1', 'C. Elegans'});
% %     set(gca, 'YScale', 'log')
% %     ylim([10e0 10e10])
%     figure;
%     [POSTHOC.COMM{i}.c, POSTHOC.COMM{i}.m, POSTHOC.COMM{i}.h, POSTHOC.COMM{i}.nms]=multcompare(ANOVA.COMM{i}.Stats,'alpha',.05/5,'ctype','bonferroni');
%     set(gca, 'XScale', 'log')
%     xlim([10e0 10e9])
%     FigHandle = setdiff(allchild(groot), FigList);
% %     xlim(
% %     if strcmp(FigHandle(j).Name,'')
% %     FBox{i}=FigHandle(2);
% %     FTable{i}=FigHandle(3);
%     FMulti{i}=FigHandle(1);
% %     end 
%     
%     
% %     set(FTable{i},'PaperPositionMode','auto');
% %     set(FTable{i},'PaperOrientation','landscape');
% %     set(FTable{i},'Position',[0 0 1920 1080]);
% %     saveas(FTable{i},['ANOVA COMM Table Across Networks sparsity ' num2str(i)],'epsc');
% %     set(FBox{i},'PaperPositionMode','auto');
% %     set(FBox{i},'PaperOrientation','landscape');
% %     set(FBox{i},'Position',[0 0 1920 1080]);
% %     saveas(FBox{i},['ANOVA COMM BoxPlot Across Networks sparsity ' num2str(i)],'epsc');
%     set(FMulti{i},'PaperPositionMode','auto');
%     set(FMulti{i},'PaperOrientation','landscape');
%     set(FMulti{i},'Position',[0 0 1920 1080]);
%     saveas(FMulti{i},['MultiComparisons COMM Plot Across Networks sparsity ' num2str(i)],'epsc');
%     clear COMMtemp 
% end
% 
%     