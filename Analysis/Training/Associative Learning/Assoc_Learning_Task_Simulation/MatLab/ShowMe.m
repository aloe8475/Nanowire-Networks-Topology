SnapInd = 10;
EdgeInd = 32;

% Uniform scales:
axesLimits.DissipationCbar = [0,5]; % (1pW*10^0 : 1pW*10^5)
axesLimits.CurrentArrowSacling = 0.12;
axesLimits.VoltageCbar = [0,Signals{1,1}(1)];
% switch Stimulus.BiasType
%     case 'DC'
%         if Stimulus.Signal(1) > 0
%             axesLimits.VoltageCbar = [0,Stimulus.Signal(1)]; % (V)
%         else
%             axesLimits.VoltageCbar = [Stimulus.Signal(1),0]; % (V)
%         end
%     case 'AC'
%         axesLimits.VoltageCbar = [min(Stimulus.Signal),max(Stimulus.Signal)]; % (V)
% end

VizSnap(snapshots{SnapInd}, Output.sources, Output.drains,Connectivity,axesLimits);
set(gcf, 'visible','on')
