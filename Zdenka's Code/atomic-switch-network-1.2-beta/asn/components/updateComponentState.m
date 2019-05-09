function updateComponentState(compPtr, dt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The function updates the state field of the input struct (which is passed
% by reference). (for example, the 'charge' fields for memristors or the
% 'filamentState' field for atomic switches.
%
% ARGUMENTS: 
% compPtr - a pointer to a struct containing the properties and current 
%           state of the electrical components of the network.
% dt - length of current timestep.
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
            compPtr.comp.charge = compPtr.comp.charge + ...
                                  (compPtr.comp.voltage./compPtr.comp.resistance) ...
                                  *dt;
        case 'atomicSwitch'
            wasOpen = abs(compPtr.comp.filamentState) >= compPtr.comp.criticalFlux;
            
            % compPtr.comp that have voltage bigger than setVoltage 
            % experience  polarity-dependent cation migration (filament can
            % either grow or shrink). The rate of change is determined by
            % the difference V-setV.
            compPtr.comp.filamentState = compPtr.comp.filamentState + ...
                                         (abs(compPtr.comp.voltage) > compPtr.comp.setVoltage) .* ...
                                         (abs(compPtr.comp.voltage) - compPtr.comp.setVoltage) .* ...
                                         sign(compPtr.comp.voltage) ...
                                         * dt;
            
            % compPtr.comp that have voltage smaller than resetVoltage
            % experience polarity-independent filament dissolution
            % (filaments always shrink). The rate of change is determined 
            % by the difference resetV-V.
            compPtr.comp.filamentState = compPtr.comp.filamentState - ...
                                         (compPtr.comp.resetVoltage > abs(compPtr.comp.voltage)) .* ...
                                         (compPtr.comp.resetVoltage - abs(compPtr.comp.voltage)) .* ...
                                         sign(compPtr.comp.filamentState) ...
                                         * dt * 10;
                                   
            % maxFlux is an upper limit on filamentState:
            compPtr.comp.filamentState (compPtr.comp.filamentState >  compPtr.comp.maxFlux) =  compPtr.comp.maxFlux(compPtr.comp.filamentState >  compPtr.comp.maxFlux);
            compPtr.comp.filamentState (compPtr.comp.filamentState < -compPtr.comp.maxFlux) = -compPtr.comp.maxFlux(compPtr.comp.filamentState < -compPtr.comp.maxFlux);
            
            % Filaments that have just disconnected suffer a blow:
            justClosed = wasOpen & (abs(compPtr.comp.filamentState) < compPtr.comp.criticalFlux);
            compPtr.comp.filamentState(justClosed) = compPtr.comp.filamentState(justClosed) / 10;
    end
end