function plotResults(resistance,current,Stimulus)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots the results of the simulation (including Fourier analysis and some
% statistics).
%
% ARGUMENTS: 
% resistance, current - the results of the simulation.
% Stimulus - the external voltage signal.
%
% OUTPUT:
% none
%
% REQUIRES:
% find_stability_histogram
% fourier_analysis
%
% USAGE:
%{
    [resistance, current, snapshots] = simulateNetwork(Equations, Components, Stimulus)
    plotResults(resistance, current, Stimulus);
%}
%
% Authors:
% Ido Marcus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cplot2 = true;
    conductance = resistance;
    figure;
    switch Stimulus.BiasType
%        case {'DC'}
        case {'DC', 'DCandWait', 'Square', 'AlonPulse'}
            subplot(1,2,1);
            % overplot C and V on same plot - use yyaxis instead of plotyy
            % (see below)
            yyaxis left; % activate left axis to plot C 
%            plot(Stimulus.TimeAxis,conductance);
            semilogy(Stimulus.TimeAxis,conductance);
            title('Conductance time series');
            %ylim([min(conductance)-.01*min(conductance), max(conductance)+0.05*max(conductance)])            
            xlim([0, max(Stimulus.TimeAxis)]);
            xlabel('Time (s)','FontSize',16);
            ylabel('Conductance (S)','FontSize',16);
            %grid on;
             yyaxis right % activate right axis to plot V
             plot(Stimulus.TimeAxis,Stimulus.Signal);
             ylim([0,3]);
             ylabel('Voltage (V)','FontSize',16);
%             % add an inset to zoom into complex voltage pulses:
%               axes('Position',[0.25,0.25,0.15,0.3]);
%               box on;
%               xlim([0,5]);
%               semilogy(Stimulus.TimeAxis,conductance);
%             plot(Stimulus.TimeAxis,Stimulus.Signal);
%             xlim([0,Stimulus.OffTime/4]);
%             ylim([0,Stimulus.AmplitudeOn]);
            
            % Fourier analysis of conductance:
            [t_freq, conductance_freq] = fourier_analysis(Stimulus.TimeAxis, conductance);
%            [t_freq, conductance_freq] = fourier_analysis(Stimulus.TimeAxis(Stimulus.TimeAxis>=1), conductance(Stimulus.TimeAxis>=1));
            % using built-in function:
%            [pwr,f] = pspectrum(conductance,Stimulus.TimeAxis,'leakage',0.5);    

            % Linear fit for log-log plot of PSD:
            fitCoef = polyfit(log10(t_freq(t_freq~=0 & t_freq<5e2)), log10(conductance_freq(t_freq~=0 & t_freq<5e2)), 1);
%            fitCoef = polyfit(log10(t_freq(t_freq>=0.2 & t_freq<=20)), log10(conductance_freq(t_freq>=0.2 & t_freq<=20)), 1);
            fitCoef(2) = 10^fitCoef(2); 
            PSDfit = fitCoef(2)*t_freq.^fitCoef(1);

            subplot(1,2,2);
            %            semilogy(t_freq,conductance_freq);
            loglog(t_freq,conductance_freq);
            xlim([min(t_freq), max(t_freq)]);
            hold on;
%            loglog(f,pwr,'g');
%            semilogy(t_freq,PSDfit,'r');
            loglog(t_freq,PSDfit,'r');
%            loglog(t_freq(t_freq>=0.2 & t_freq<=20),PSDfit(t_freq>=0.2 & t_freq<=20),'r');
            text(0.5,0.8,sprintf('\\beta=%.1f', -fitCoef(1)),'Units','normalized','Color','r','FontSize',18);
            title(strcat('Conductance PSD '));
            xlabel('Frequency (Hz)');
            ylabel('PSD');
            ylim([min(conductance_freq)/10,max(conductance_freq)*10]);
            set(gca,'Ytick',10.^(-20:1:20));
            %set(gca,'Xtick',10.^(-20:1:20));
            grid on;
            
            % Plot spectrogram:
%            figure('Name','Spectrogram');
%            pspectrum(conductance,Stimulus.TimeAxis,'leakage',0.5,'TimeResolution',0.5,'OverlapPercent',60,'spectrogram');
            
%             subplot(2,2,3);
            % PDF of duration of time interval between changes:
%             [stability_hist_bin, stability_hist_values] = find_stability_histogram(conductance,Stimulus.dt);
%             
%             % Linear fit of this PDF in log-log plot:
%             fitCoef2 = polyfit(log10(stability_hist_bin(stability_hist_bin~=0 & stability_hist_values~=0)), log10(stability_hist_values(stability_hist_bin~=0 & stability_hist_values~=0)), 1);
%             fitCoef2(2) = 10^fitCoef2(2); 
%             PDFfit = fitCoef2(2)*stability_hist_bin.^fitCoef2(1);
%             plot(stability_hist_bin,stability_hist_values,'*');
%             title('PDF of length of interval between changes');
%             xlabel('Time interval (sec)')
%             ylabel('PDF');
%             grid on;
% 
%             subplot(2,2,4);
%             loglog(stability_hist_bin,stability_hist_values,'*');
%             hold on
%             loglog(stability_hist_bin,PDFfit,'r');
%             text(0.5,0.8,sprintf('\\beta=%.1f', -fitCoef2(1)),'Units','normalized','Color','r','FontSize',18);
%             title('PDF of length of interval between changes (log-log)');
%             xlabel('Time interval (sec)')
%             ylabel('PDF');
%             xlim([min(stability_hist_bin),max(stability_hist_bin)]);
%             grid on;
            
        case {'AC', 'ACsaw'}
            subplot(2,1,1);
            [t_freq, extCur_freq] = fourier_analysis(Stimulus.TimeAxis, current);
            % use yyaxis instead of plotyy
            %[axesHandles, curLine, volLine] = plotyy(Stimulus.TimeAxis*Stimulus.Frequency,current,...
                                                     %Stimulus.TimeAxis*Stimulus.Frequency,Stimulus.Signal);
            yyaxis left % activate left axis to plot V
            plot(Stimulus.TimeAxis,Stimulus.Signal);
            ylim([-5,5]);
            ylabel('Voltage (V)','FontSize',14);
            yyaxis right; % activate right axis to plot I 
            plot(Stimulus.TimeAxis,current);
            %ylim([-2e-7,2e-7]);
            ylabel('Current (A)','FontSize',14);
%            ylabel(axesHandles(1),'Voltage (V)');
%            ylabel(axesHandles(2),'Current (A)');
            xlabel('Cycle Number','FontSize',14);
            title('Current vs voltage sweep cycle');
%             curLine.Color='red';
%             curLine.LineWidth = 1;
%             volLine.Color='cyan';
%             volLine.LineWidth = 0.1;
            grid on;

            subplot(2,1,2);
            
%             hermonicsVecor = t_freq/Stimulus.Frequency;
%             extCur_freq = extCur_freq(hermonicsVecor<10);
%             hermonicsVecor = hermonicsVecor(hermonicsVecor<10);
%             plot(hermonicsVecor,extCur_freq);
%             xlabel('Harmonic')
%             ylabel('Power density')
%             title('Power density spectrum of the current')
%             grid on;
            
            % Plot I-V curve:
%            figure('Name','I-V curve');
            plot(Stimulus.Signal,current);
            xlabel('V','FontSize',14);
            ylabel('I','FontSize',14);
            
         case {'Ramp'}
             figure
%              % plot C and V in separate plots:
%              subplot(2,1,1);
%              plot(Stimulus.TimeAxis,conductance);
%              title('Network conductance time series')
%              xlabel('Time (sec)');
%              subplot(2,1,2);
%              plot(Stimulus.TimeAxis,Stimulus.Signal);
%              title('Input voltage signal')
%              xlabel('Time (sec)');
%              ylabel('Voltage (V)');
             % plot C vs V and I vs V in same plots:
             subplot(2,1,1);
             yyaxis left; % activate left axis to plot C 
             semilogy(Stimulus.TimeAxis,conductance);
             title('Conductance time series');
             xlim([0, max(Stimulus.TimeAxis)]);
             xlabel('Time (sec)','FontSize',16);
             ylabel('Conductance (S)','FontSize',16);
             %grid on;
             yyaxis right % activate right axis to plot V
             plot(Stimulus.TimeAxis,Stimulus.Signal);
             ylim([0,5]);
             ylabel('Voltage (V)','FontSize',16);
             subplot(2,1,2);
             yyaxis left; % activate left axis to plot I 
             semilogy(Stimulus.TimeAxis,current);
             title('Current time series');
             xlim([0, max(Stimulus.TimeAxis)]);
             xlabel('Time (sec)','FontSize',16);
             ylabel('Current (A)','FontSize',16);
             %grid on;
             yyaxis right % activate right axis to plot V
             plot(Stimulus.TimeAxis,Stimulus.Signal);
             ylim([0,5]);
             ylabel('Voltage (V)','FontSize',16);
            % add an inset to zoom into complex voltage pulses:
%               axes('Position',[0.58,0.6,0.3,0.3]);
%               box on;
%               %xlim([0,5]);
%               semilogy(Stimulus.TimeAxis,conductance);
%             plot(Stimulus.TimeAxis,Stimulus.Signal);
%             xlim([0,Stimulus.OffTime/4]);
%             ylim([0,Stimulus.AmplitudeOn]);

        case {'Ramp', 'AC', 'ACsaw'}
            figure
            [axesHandles, ~, ~] = plotyy(Stimulus.TimeAxis,conductance,...
                                         Stimulus.TimeAxis,Stimulus.Signal);
            title('External voltage and conductance as function of time')
            xlabel('Time (sec)');
            ylabel(axesHandles(1),'Conductance (S)');
            ylabel(axesHandles(2),'External Voltage (V)');
            
            
        case 'Triangle'
            subplot(121)
            [axesHandles, ~, ~] = plotyy(Stimulus.TimeAxis,conductance,...
                                         Stimulus.TimeAxis,Stimulus.Signal);
            title('External voltage and conductance as function of time')
            xlabel('Time (sec)');
            ylabel(axesHandles(1),'Conductance (S)');
            ylabel(axesHandles(2),'External Voltage (V)');
            subplot(122)
            plot(Stimulus.Signal, current)
            xlabel('Bias [V]');
            ylabel('Current [A]');
            title('Applied Voltage vs Output Current') 
    end
    
    %{
    if cplot2
        figure
    %              % plot C and V in separate plots:
    %              subplot(2,1,1);
    %              plot(Stimulus.TimeAxis,conductance);
    %              title('Network conductance time series')
    %              xlabel('Time (sec)');
    %              subplot(2,1,2);
    %              plot(Stimulus.TimeAxis,Stimulus.Signal);
    %              title('Input voltage signal')
    %              xlabel('Time (sec)');
    %              ylabel('Voltage (V)');
         % plot C vs V and I vs V in same plots:
         subplot(2,1,1);
         yyaxis left; % activate left axis to plot C 
         semilogy(Stimulus.TimeAxis,conductance);
         title('Conductance time series');
         xlim([0, max(Stimulus.TimeAxis)]);
         xlabel('Time (sec)','FontSize',16);
         ylabel('Conductance (S)','FontSize',16);
         %grid on;
         yyaxis right % activate right axis to plot V
         plot(Stimulus.TimeAxis,Stimulus.Signal);
         ylim([0,5]);
         ylabel('Voltage (V)','FontSize',16);
         subplot(2,1,2);
         yyaxis left; % activate left axis to plot I 
         semilogy(Stimulus.TimeAxis,current);
         title('Current time series');
         xlim([0, max(Stimulus.TimeAxis)]);
         xlabel('Time (sec)','FontSize',16);
         ylabel('Current (A)','FontSize',16);
         %grid on;
         yyaxis right % activate right axis to plot V
         plot(Stimulus.TimeAxis,Stimulus.Signal);
         ylim([0,5]);
         ylabel('Voltage (V)','FontSize',16);
        % add an inset to zoom into complex voltage pulses:
    %               axes('Position',[0.58,0.6,0.3,0.3]);
    %               box on;
    %               %xlim([0,5]);
    %               semilogy(Stimulus.TimeAxis,conductance);
    %             plot(Stimulus.TimeAxis,Stimulus.Signal);
    %             xlim([0,Stimulus.OffTime/4]);
    %             ylim([0,Stimulus.AmplitudeOn]);
    end
    %}
    
end