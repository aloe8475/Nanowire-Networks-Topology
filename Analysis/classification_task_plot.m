%% CLASSIFICATION ANALYSIS
close all

type='d'; %s - same, d - different

if type =='d'
    filepath='C:\Users\aloe8475\Dropbox (Sydney Uni)\Data\ASN_simulation\Python\ASN\data\Collapsed\Variable Path Length\';
        filename='exp4_classification_v1_variablePathLength_sp1000.mat';
     save_directory='C:\Users\aloe8475\Documents\PhD\GitHub\CODE\Data\Figures\Explore Analysis\Classification Analysis\Variable Path Lengths\';

elseif type =='s'
    filepath= 'C:\Users\aloe8475\Dropbox (Sydney Uni)\Data\ASN_simulation\Python\ASN\data\Collapsed\Path Length 4\';
    filename='exp4_classification_v1_PathLength4_sp1000.mat';
    save_directory='C:\Users\aloe8475\Documents\PhD\GitHub\CODE\Data\Figures\Explore Analysis\Classification Analysis\Path Length 4\';

end
load([filepath filename]);

%Reshape data
for i = 1:length(Data)
    test.accuracy(i,:)=[Data{i}.test1_accuracy];
    train.accuracy(i,:)=[Data{i}.train_accuracy];
    test.class(i,:)=[Data{i}.test1_class];
    train.class(i,:)=[Data{i}.train_class];
    test.sep(i,:)=[Data{i}.test1_sep];
train.sep(i,:)=[Data{i}.train_sep];
    test.pathLength(i,:)=[Data{i}.test1_pathLength];
    train.pathLength(i,:)=[Data{i}.train_pathLength];
end 
%% ANOVAS

        FigList = allchild(groot);
        
[ANOVA.p,ANOVA.AnovaTab,ANOVA.Stats] = anovan(test.accuracy,{train.class,test.class,train.sep,test.sep},'varnames',strvcat('Training Class', 'Testing Class', 'Training TimeDelay', 'Testing TimeDelay'));
%     [ANOVA.p,ANOVA.AnovaTab,ANOVA.Stats] = anovan(exp4classification.test1_accuracy,{exp4classification.train_sep,exp4classification.test1_sep,exp4classification.train_pathLength,exp4classification.test1_pathLength},'varnames',strvcat('Training TimeDelay', 'Testing TimeDelay','Training PathLength','Testing PathLength'));
figure;
[POSTHOC.c, POSTHOC.m, POSTHOC.h, POSTHOC.nms] = multcompare(ANOVA.Stats,'Dimension',[3,4],'ctype','bonferroni'); %dimension 1 = training class, 2 = testing class, 3 = training timedelay, 4 = testing timedelay 
   
FigHandle = setdiff(allchild(groot), FigList);
    Ftable=FigHandle(1);
    FMulti=FigHandle(2);
    
    clear FigHandle FigList

%SAVE FIGURES:

saveas(Ftable,[save_directory '100nw_LDA_Classification_ANOVA_Table' date],'png');
print(Ftable,'-painters','-dpdf','-bestfit','-r300',[save_directory '100nw_LDA_Classification_ANOVA_Table_' date '.pdf']);
set(FMulti,'PaperPositionMode','auto');
set(FMulti,'PaperOrientation','landscape');
set(FMulti,'Position',[0 0 1920 1080]);
saveas(FMulti,[save_directory '100nw_LDA_Classification_MultipleComparisons_' date],'png');
print(FMulti,'-painters','-dpdf','-bestfit','-r300',[save_directory '100nw_LDA_Classification_MultipleComparisons_' date '.pdf']);
   