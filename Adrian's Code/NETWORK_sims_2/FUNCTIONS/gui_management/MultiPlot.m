function handles=MultiPlot(handles,IndexNet,IndexSim,TimeStep)
%% The most important function for gui management
% Takes the selected network and simulation and makes a multiplot
% in the main axes of the figure, with the desired kind of plot
% selected in the PlotList
% there are 4 different types of plot> timeseries,network, frequency, graph
% One a subplot is created, the for axes are stored in the main_gui handles
% struct (output function of multiplot) as a handles.multiAx cell array.
% the handles.Status property controls wether network is redrawn or only
% updated relevant information (such as new curents or voltages) are drawn.


%% relevant variables from handles struct


Network=handles.Networks{IndexNet};
Sim=Network.Simulations{IndexSim};

PlotSels=handles.PlotSel;
Status=handles.PlotStatus; %first,redraw

NoN=(strcmp(PlotSels.Type,'NoPlot'));
EffPlot=1:4;
EffPlot=EffPlot(~NoN);
[OuterPosArr,dim]=SetPosSubPlots(sum(NoN));
if isequal(Status,'first') handles.multiAx={}; end



for i=1:length(EffPlot)
    if isequal(Status,'first')
        handles.multiAx{i}=subplot(dim(1),dim(2),i);
        currAx=gca;
        set(currAx,'OuterPosition',OuterPosArr(i,:));
        switch PlotSels.Name{EffPlot(i)}
            case {'TimeSeries','FrequencySeries'}
            otherwise
                % matlab subplot sets a very wide margin between plots, this
                % function shrinks this margin to make the figures wider.
                set(currAx,'Position',getPosTight(OuterPosArr(i,:),currAx.TightInset));
                
        end
        
    else
        currAx=handles.multiAx{i};
    end
    
    funcHand=str2func(strcat('Plot',PlotSels.Name{EffPlot(i)}));
    funcHand(currAx,handles.NodeList,Sim,PlotSels.Type{EffPlot(i)},Status,TimeStep,...
        handles.GraphPop.String{handles.GraphPop.Value});
    UpdatePlotSettings(currAx,PlotSels.Name{EffPlot(i)},PlotSels.Type{EffPlot(i)});
    
end



end