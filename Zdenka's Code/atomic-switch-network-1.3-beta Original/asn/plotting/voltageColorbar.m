function voltageColorbar(biasMode)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script plots a voltage colorbar (green to red), that is to be 
% manually pasted onto a snapshot figure.
%
% ARGUMENTS: 
% biasMode - 'DC' or 'AC'.
%
% OUTPUT:
% none.
%
% REQUIRES:
% none.
%
% USAGE:
%{
    voltageColorbar('DC');
%}
%
% Author:
% Ido Marcus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    close all;
    figure;
    colormap([(1:-0.01:0)', (0:0.01:1)', zeros(length(0:0.01:1),1)]);
    vColorbar = colorbar('southoutside');
    switch biasMode
        case 'DC'
            caxis([0,1]);
            vColorbar.Ticks = [0,1];
            vColorbar.TickLabels = {'0', 'V_e_x_t'};
        case 'AC'
            caxis([-1,1]);
            vColorbar.Ticks = [-1,0,1];
            vColorbar.TickLabels = {'-V_e_x_t', '0', 'V_e_x_t'};
    end
    vColorbar.Label.String = 'Voltage';
    set(gca,'FontSize',50);
    axis off;
end