function snapshotFigure = VizSnap(snapshot, sources, drains, connectivity, axesLimits)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This funciton generates a visualization of a snapshot of the network.
% This includes: the spatial distribution of wires, the location of
% contacts, the state of each switch (ON\OFF), the voltage on each wire the
% currents distribution and the dissipation (heat-map).
%
% ARGUMENTS: 
% snapshot - A struct containing the voltage, resistance, etc. of the 
%            electrical components in the network, at a particular 
%            timestamp.
% contact - the indices of the two wires that serve as contacts.
% connectivity - A structure with the adjacency matrix as well as the
%                spatial information of the wires (location, orientation
%                etc.).
% whatToPlot - A structure of boolean flags, controlling the contents of
%              the output figure. Fields:
%                .Nanowires
%                .Contacts
%                .Dissipation
%                .Currents
%                .Voltages
%
% axesLimits - A structure of limits for the different axes used. This is
%              usefull in order to compile a movie, in which the colorbar,
%              length of arrows etc. should be consistent between the
%              different frames. Fields:
%                .DissipationCbar - since the colorbar is shown in
%                                   logarithmic scale, its limits will be:
%                                   [10^axesLimits.DissipationCbar(1),
%                                   10^axesLimits.DissipationCbar(2)]     
%                                   (in pW).
%                .CurrentArrowSacling - a common factor that scales all the
%                                       current arrows.
%                .VoltageCbar - limits of the voltage colorbar. for
%                               (positive DC bias that should be [0,Vext].
%
% OUTPUT:
% snapshotFigure - a figure object containing the wanted map of the
%                  network.
%
% REQUIRES:
% getAbsoluteVoltage
%
% Authors:
% Ido Marcus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    contacts = [sources, drains];
    whatToPlot = struct(...
                    'Nanowires',    true, ...
                    'Contacts',     true, ...
                    'Dissipation',  true, ...
                    'Currents',     false, ...
                    'Voltages',     true  ...
                    );



    %% input control
    if ~strcmp(connectivity.WhichMatrix, 'nanoWires')
        error('Cannot visualize a snapshot for connectivity graphs that have no spatial meaning');
    end

    %% preliminaries
    snapshotFigure = figure('visible','off','Position', [0, 0, 1000, 550]);
    set(gca,'Color',[0.2 0.2 0.2]);
    hold all;
        
    %% nanowires and voltage distribution:  
    if whatToPlot.Voltages    
        % compute absolute voltages:
        absoluteVoltage = getAbsoluteVoltage(snapshot, connectivity, contacts);
        
        % trim results so that they don't saturate the color-scale:
        absoluteVoltage(absoluteVoltage < axesLimits.VoltageCbar(1)) = axesLimits.VoltageCbar(1);
        absoluteVoltage(absoluteVoltage > axesLimits.VoltageCbar(2)) = axesLimits.VoltageCbar(2);
        
        % linearly transform results to the range [0,1]:
        voltageColorCode = (absoluteVoltage - axesLimits.VoltageCbar(1)) / (axesLimits.VoltageCbar(2) - axesLimits.VoltageCbar(1));
        
        % construct RGB triplets (from green near contact(1) to red near contact(2)):
        voltageColorCode = [1-voltageColorCode,voltageColorCode,zeros(connectivity.NumberOfNodes,1)];
        
        % add color bar:
            % Matlab allows only one colorbar per axes, see
            % voltageColorbar.
    end
    
    if whatToPlot.Nanowires
        if whatToPlot.Voltages    
            % nanowires with a voltage color-code
            lineColor = voltageColorCode;
        else
            % only (white) nanowires
            lineColor = ones(connectivity.NumberOfNodes,3);
        end
        for currWire=1:connectivity.NumberOfNodes
                line([connectivity.WireEnds(currWire,1),connectivity.WireEnds(currWire,3)],[connectivity.WireEnds(currWire,2),connectivity.WireEnds(currWire,4)],'Color',lineColor(currWire,:),'LineWidth',0.5)
        end
    else
        if whatToPlot.Voltages
            % no nanowires, only points at the centers of the nanowires with a voltage color-code
            markerSize = 100;
            scatter(connectivity.VertexPosition(:,1), ...
                    connectivity.VertexPosition(:,2), ...
                    markerSize,                       ...
                    voltageColorCode,                 ...
                    'filled',                         ...
                    's'                               ...
                    );
                    
        else
            % nothing.
        end
    end
    
    %% dissipation:
    if whatToPlot.Dissipation
        junctionSize = 50;
        
        % calculate power consumption:
        power = 1e9*(snapshot.Voltage).^2./(snapshot.Resistance); % Joule heating is V^2/R.
        
        % if possible, give different marker to open switches, otherwise
        % just plot all the switches in the same manner:
        if isfield(snapshot, 'OnOrOff')
            on = snapshot.OnOrOff(1:end-1);

            onDots = scatter(connectivity.EdgePosition(on,1),    ...
                             connectivity.EdgePosition(on,2),    ...
                             1.5*junctionSize,                       ...
                             log10(power(on)),                     ...
                             'filled',                           ...
                             'o');

            offDots = scatter(connectivity.EdgePosition(~on,1),   ...
                              connectivity.EdgePosition(~on,2),   ...
                              junctionSize,                       ...
                              log10(power(~on)),                    ...
                              'filled',                           ...
                              's');

            % legend:
            if any(on) && any(~on)
                leg = legend([onDots,offDots],{'ON switch','OFF switch'});
            elseif any(on)
                leg = legend(onDots,{'ON switch'});
            else
                leg = legend(offDots,{'OFF switch'});
            end
            leg.Color = 'white';

        else
            scatter(connectivity.EdgePosition(:,1),     ...
                    connectivity.EdgePosition(:,2),     ...
                    junctionSize,                       ...
                    log10(power),                         ...
                    'filled');
        end
        
        % colorbar:               
        dissipationColorbar  = colorbar;
        dissipationColorbar.Label.String = 'Power (logarithmic scale)';
        %upperLimit = ceil(max(log10(power)));
        upperLimit = axesLimits.DissipationCbar(2);
        caxis([0,upperLimit]);
        dissipationColorbar.Ticks = 0:1:upperLimit;
        dissipationColorbar.TickLabels(:) = '';
            % tick labels:
            labelBank = {'pW', '', '', 'nW', '', '', '\muW', '', '', 'mW', '', '', 'W', '', '', 'kW', '', '', 'MW'};
            dissipationColorbar.TickLabels = labelBank(1:upperLimit+1);            
        colormap('parula');
    end
    
    %% currents:
    if whatToPlot.Currents
        % Allocate space (assuming no intersections at ends of wires):
        numberOfSections = 2*length(connectivity.EdgeList) - connectivity.NumberOfNodes + 2;
        sectionCenterX  = zeros(numberOfSections,1);
        sectionCenterY  = zeros(numberOfSections,1);
        sectionCurrentX = zeros(numberOfSections,1);
        sectionCurrentY = zeros(numberOfSections,1);
        numSectionsDone = 0;
        
        % Calculate currents:
        currents = 1e6*snapshot.Voltage./snapshot.Resistance; % (nA)
        
        % Calculate wire angles ([-pi/2,pi/2]):
                % first [0,pi]
        wireAngles = mod(atan2(connectivity.WireEnds(:,4)-connectivity.WireEnds(:,2), connectivity.WireEnds(:,3)-connectivity.WireEnds(:,1)),pi);
            % The modulu operation makes sure that the result is between
            % 0 and pi. It is just for safety, since the ends are sorted
            % such that WireEnds(:,4)>WireEnds(:,2).        
                % then [-pi/2,pi/2]
        wireAngles(wireAngles>pi/2) = wireAngles(wireAngles>pi/2) - pi;
            % Since the sections will soon be sorted from left to right, a
            % positive current value along a section must always have a
            % positive cosine value, and vise versa.
   
        % Wire by wire:
        sourcePoints = [];
        drainPoints = [];
        for currWire=1 : connectivity.NumberOfNodes        
 
            % Find the indices of edges (=intersections) relevant for this
            % vertex (=wire):
            relevantEdges = find(connectivity.EdgeList(1,:) == currWire | connectivity.EdgeList(2,:) == currWire);
            
            % Sort them according to physical location (left-to-right, or
            % if the wire is vertical then bottom-up):
            if connectivity.WireEnds(currWire,1) ~= connectivity.WireEnds(currWire,3)
                [~,I] = sort(connectivity.EdgePosition(relevantEdges,1));
            else
                [~,I] = sort(connectivity.EdgePosition(relevantEdges,2));
            end
            relevantEdges = relevantEdges(I);
            
            % Calculate the current along each section of the wire:
            direction = ((currWire ~= connectivity.EdgeList(1,relevantEdges))-0.5)*2;       
                % Using the convention that currents are defined to flow
                % from wires with lower index to wires with higher index, 
                % and that in the field EdgeList the upper row always 
                % contains lower indices. 
            wireCurrents = cumsum(currents(relevantEdges(1:end-1)).*direction(1:end-1)'); 
                % The first element in wireCurrents is the current in the
                % section between relevantEdge(1) and relevantEdge(2). We
                % assume that between relevantEdge(1) and the closest wire
                % end there's no current. Then, acording to a KCL
                % equation, there's also no current from relevantEdge(end)
                % to the other wire end (that's why the end-1 in the 
                % expression). 
                % The only two exceptions are the two contacts, where the 
                % contact point is defined as the wire end closest to 
                % relevantEdge(end).
            
            % Accumulate for a quiver (vector field) plot:
            first = numSectionsDone + 1;
            last = first + length(wireCurrents) - 1;
            sectionCenterX(first:last)  = mean([connectivity.EdgePosition(relevantEdges(1:end-1),1), connectivity.EdgePosition(relevantEdges(2:end),1)],2);
            sectionCenterY(first:last)  = mean([connectivity.EdgePosition(relevantEdges(1:end-1),2), connectivity.EdgePosition(relevantEdges(2:end),2)],2);
            sectionCurrentX(first:last) = cos(wireAngles(currWire))*wireCurrents;
            sectionCurrentY(first:last) = sin(wireAngles(currWire))*wireCurrents;
            numSectionsDone = last;
            
            % Contacts stuff:
            if any([sources, drains] == currWire)
                % Find the relevant end of the wire:
                if connectivity.WireEnds(currWire,1) ~= connectivity.WireEnds(currWire,3)
                    if connectivity.WireEnds(currWire,1) < connectivity.WireEnds(currWire,3)
                        contactEnd = connectivity.WireEnds(currWire,3:4);
                    else
                        contactEnd = connectivity.WireEnds(currWire,1:2);
                    end
                else
                    if connectivity.WireEnds(currWire,2) < connectivity.WireEnds(currWire,4)
                        contactEnd = connectivity.WireEnds(currWire,3:4);
                    else
                        contactEnd = connectivity.WireEnds(currWire,1:2);
                    end
                end

                % Add total current arrow:
                totalCurrent = sum(currents(relevantEdges).*direction');
                sectionCenterX(numSectionsDone+1)  = mean([connectivity.EdgePosition(relevantEdges(end),1), contactEnd(1)],2);
                sectionCenterY(numSectionsDone+1)  = mean([connectivity.EdgePosition(relevantEdges(end),2), contactEnd(2)],2);
                sectionCurrentX(numSectionsDone+1) = cos(wireAngles(currWire))*totalCurrent; 
                sectionCurrentY(numSectionsDone+1) = sin(wireAngles(currWire))*totalCurrent;
                numSectionsDone = numSectionsDone + 1;
                
                % Gather contact point location, if needed:
                if whatToPlot.Contacts
                % The location of the current wire's end which is 
                % rightmost (or if vertical, upper) is defined as the 
                % contact POINT:
                    if any(currWire == sources)
                        sourcePoints = [sourcePoints;contactEnd];
                    end
                    if any(currWire == drains)
                        drainPoints  = [drainPoints;contactEnd];
                    end
                end
            end
        end
        
        % Plot current arrows:
        quiver(sectionCenterX,sectionCenterY,sectionCurrentX/axesLimits.CurrentArrowSacling,sectionCurrentY/axesLimits.CurrentArrowSacling,0,'Color','r','LineWidth',1.5);
        %quiver(sectionCenterX,sectionCenterY,sectionCurrentX,sectionCurrentY,'Color','w','LineWidth',1);
    end
    
    %% contacts:   
    if whatToPlot.Contacts
        if whatToPlot.Currents
            scatter(sourcePoints(:,1),sourcePoints(:,2),200,[0 1 0],'filled','h');
            scatter(drainPoints(:,1),drainPoints(:,2),200,[1 0 0],'filled','h');
            sourceTextPosition = sourcePoints;
            drainTextPosition  = drainPoints;
        else
            line([connectivity.WireEnds(contacts(1),1),connectivity.WireEnds(contacts(1),3)],[connectivity.WireEnds(contacts(1),2),connectivity.WireEnds(contacts(1),4)],'Color','g','LineWidth',0.2)
            line([connectivity.WireEnds(contacts(2),1),connectivity.WireEnds(contacts(2),3)],[connectivity.WireEnds(contacts(2),2),connectivity.WireEnds(contacts(2),4)],'Color','r','LineWidth',0.2)
            sourceTextPosition = connectivity.VertexPosition(contacts(1),:);
            drainTextPosition  = connectivity.VertexPosition(contacts(2),:);
        end
        text(sourceTextPosition(:,1), sourceTextPosition(:,2), '  (S)',  'Color', 'g', 'FontSize', 16);
        text(drainTextPosition(:,1),  drainTextPosition(:,2),  '  (D)',  'Color', 'r', 'FontSize', 16);
    end
    
%      scatter(connectivity.EdgePosition(EdgeInd,1),    ...
%              connectivity.EdgePosition(EdgeInd,2),    ...
%              75,                    ...        
%              [1 0.65 0],                            ...
%              'filled',                           ...
%              '^');
    
    %% title, axes labels and limits:
    title(sprintf('t=%.2f (sec)', snapshot.Timestamp));
    xlabel('x (\mum)');
    ylabel('y (\mum)');
    shoulder = 800;
    axis([-shoulder,connectivity.GridSize(1)+shoulder,-shoulder,connectivity.GridSize(2)+shoulder]);
    axis square;
end