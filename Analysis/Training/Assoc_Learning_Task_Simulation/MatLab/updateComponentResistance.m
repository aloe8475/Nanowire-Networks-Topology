function resistance = updateComponentResistanceMulti(compPtr)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The function updates the 'resistance' field of the input struct (which is 
% passed by reference).
%
% ARGUMENTS: 
% compPtr - a pointer to a struct containing the properties and current 
%           state of the electrical components of the network.
%
% OUTPUT:
% none
%
% REQUIRES:
% none
%
% Authors:
% Ido Marcus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    switch compPtr.comp.type
        case 'memristor'
            % Relevant fields: 
                % charge - time integral over current for every component.
                % memristance properties - details of the (... _/-\_/-\ ...) shape.
                  
            % Even function of charge:
            charge = abs(compPtr.comp.charge);
            
            % Periodic function of charge:
            %charge = bsxfun(@mod, charge, Components.period);

            % Create a _/-\_ shaped memristance (each component with its own
            % specific values):
            resistance = (charge <= compPtr.comp.lowThreshold).*compPtr.comp.offResistance;

            resistance = ...
                resistance + (charge > compPtr.comp.lowThreshold & ...
                              charge <= compPtr.comp.highThreshold).* ...
                (  ...
                   compPtr.comp.offResistance + ...
                   ((compPtr.comp.onResistance-compPtr.comp.offResistance)./(compPtr.comp.highThreshold-compPtr.comp.lowThreshold)) .* ...
                   (charge-compPtr.comp.lowThreshold) ...
                );

            resistance = resistance + (charge >  compPtr.comp.highThreshold).*compPtr.comp.onResistance;

        case 'atomicSwitch'
            compPtr.comp.OnOrOff = abs(compPtr.comp.filamentState) >= compPtr.comp.criticalFlux;
            compPtr.comp.OnOrOff(compPtr.comp.identity == 0) = true;        
            % passive elements (resistors) are allways considered as "open" switches
            
            resistance = (~compPtr.comp.OnOrOff) .* compPtr.comp.offResistance;
            resistance = resistance + ...
                         ( compPtr.comp.OnOrOff) .* compPtr.comp.onResistance;
        case 'resistor'
            resistance = zeros(size(compPtr.comp.identity)); 
            % That's a place-holder. If resistance is not initialized, the 
            % next statement which takes care of passive elements creates a
            % row rather than a column vector.
    end    
    
    % Components that are resistors have resistance 'onResistance',
    % regardless of anything else:
    resistance(compPtr.comp.identity == 0) = compPtr.comp.onResistance(compPtr.comp.identity == 0);
    
    % Modify the input with the updated resistance values:
    compPtr.comp.resistance = resistance;
end