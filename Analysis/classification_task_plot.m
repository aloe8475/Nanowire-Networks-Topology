%% CLASSIFICATION ANALYSIS
filename='exp4_classification_v1.5_PathLength4_sp1000.mat';
filepath= 'C:\Users\aloe8475\Dropbox (Sydney Uni)\Data\ASN_simulation\Python\ASN\data\Collapsed\Path Length 4\';
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
[ANOVA.p,ANOVA.AnovaTab,ANOVA.Stats] = anovan(test.accuracy,{train.class,test.class,train.sep,test.sep},'varnames',strvcat('Training Class', 'Testing Class', 'Training TimeDelay', 'Testing TimeDelay'));
%     [ANOVA.p,ANOVA.AnovaTab,ANOVA.Stats] = anovan(exp4classification.test1_accuracy,{exp4classification.train_sep,exp4classification.test1_sep,exp4classification.train_pathLength,exp4classification.test1_pathLength},'varnames',strvcat('Training TimeDelay', 'Testing TimeDelay','Training PathLength','Testing PathLength'));
mutipleComparisons = multcompare(ANOVA.Stats,'Dimension',[1,3,4],'ctype','bonferroni'); %dimension 1 = training class, 2 = testing class, 3 = training timedelay, 4 = testing timedelay