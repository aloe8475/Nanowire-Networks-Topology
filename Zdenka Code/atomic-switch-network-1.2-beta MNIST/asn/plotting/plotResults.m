function plotResults(resistance,current,Stimulus)
global plotpsd;
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
    
    conductance = 1./resistance;
    
    switch Stimulus.BiasType
        case {'DC', 'DCandWait'}
            figure;
%           subplot(1,2,1);

            % overplot C and V on same plot - use yyaxis instead of plotyy

            % (see below)

            yyaxis left; % activate left axis to plot C 

            plot(Stimulus.TimeAxis,conductance);

%            semilogy(Stimulus.TimeAxis,conductance);

            title('Conductance time series');

            ylim([min(conductance)-.01*min(conductance), max(conductance)+0.05*max(conductance)])            

            xlim([0, max(Stimulus.TimeAxis)]);

            xlabel('Time (s)','FontSize',16);
%             set(gca,'XTick',[0:10:120]);
            set(gca,'XMinorTick','on','YMinorTick','on');
%             set(gca,'XTickLabel',{'0','','10','','20','','30','','40','','50','','60','','70','','80','','90','','100','','110','','120'});
            ylabel('Conductance (S)','FontSize',16);

            %grid on;

             yyaxis right % activate right axis to plot V

             plot(Stimulus.TimeAxis,Stimulus.Signal);

             ylim([0,3]);

             ylabel('Voltage (V)','FontSize',16);
             set(gca,'FontSize',16);  set(gca,'XMinorTick','on','YMinorTick','on');set(gca,'TickLength',[0.02, 0.01]);
                
            % Fourier analysis of conductance:
            [t_freq, conductance_freq] = fourier_analysis(Stimulus.TimeAxis, conductance);
            
            % Linear fit for log-log plot of PSD:
            fitCoef = polyfit(log10(t_freq(t_freq~=0 & t_freq<1e2)), log10(conductance_freq(t_freq~=0 & t_freq<1e2)), 1);
            fitCoef(2) = 10^fitCoef(2); 
            PSDfit = fitCoef(2)*t_freq.^fitCoef(1);

            % PDF of duration of time interval between changes:
            [stability_hist_bin, stability_hist_values] = find_stability_histogram(conductance,Stimulus.dt);
            
            % Linear fit of this PDF in log-log plot:
            fitCoef2 = polyfit(log10(stability_hist_bin(stability_hist_bin~=0 & stability_hist_values~=0)), log10(stability_hist_values(stability_hist_bin~=0 & stability_hist_values~=0)), 1);
            fitCoef2(2) = 10^fitCoef2(2); 
            PDFfit = fitCoef2(2)*stability_hist_bin.^fitCoef2(1);
            
% %             subplot(1,2,1);
%             plot(Stimulus.TimeAxis,conductance);
%             title('Conductance as function of time');
%             ylim([min(conductance)-.01*min(conductance), max(conductance)+0.05*max(conductance)])
%             xlabel('Time (sec)');
%             ylabel('Conductance (S)');
%             grid on;

%             subplot(1,2,2);
global plotpsd;
if plotpsd
              figure;
            semilogy(t_freq,conductance_freq);
            hold on;
            semilogy(t_freq,PSDfit,'r');
            text(0.5,0.8,sprintf('\\beta=%.1f', -fitCoef(1)),'Units','normalized','Color','r','FontSize',18);
            title('Power density spectrum of the conductance');
            xlabel('Frequency (Hz)');
            ylabel('Power density');
            xlim([0, 100]);
            ylim([min(conductance_freq)/10,max(conductance_freq)*10]);
            set(gca,'Ytick',10.^(-20:1:20));
            %set(gca,'Xtick',10.^(-20:1:20));
            set(gca,'FontSize',16);  set(gca,'XMinorTick','on','YMinorTick','on');
            grid on;
end

%             subplot(2,2,3);
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
            
        case 'AC'
%             [t_freq, extCur_freq] = fourier_analysis(Stimulus.TimeAxis, conductance);
            figure;
%           subplot(1,2,1);

            % overplot C and V on same plot - use yyaxis instead of plotyy

            % (see below)

            yyaxis left; % activate left axis to plot C 

            plot(Stimulus.TimeAxis,current(:,1));

%            semilogy(Stimulus.TimeAxis,conductance);

            title('Conductance time series');

%             ylim([min(conductance)-.01*min(conductance), max(conductance)+0.05*max(conductance)])            

            xlim([0, max(Stimulus.TimeAxis)]);

            xlabel('Time (s)','FontSize',16);

            ylabel('Conductance (S)','FontSize',16);
           set(gca,'FontSize',16);  set(gca,'XMinorTick','on','YMinorTick','on');
            %grid on;

             yyaxis right % activate right axis to plot V

             plot(Stimulus.TimeAxis,Stimulus.Signal);

             ylim([-3,3]);

             ylabel('Voltage (V)','FontSize',16);
            set(gca,'FontSize',16);  set(gca,'XMinorTick','on','YMinorTick','on');
%             subplot(2,1,1);
if 0
 figure;
%             [axesHandles, curLine, volLine] = plotyy(Stimulus.TimeAxis*Stimulus.Frequency,current,...
%                                                      Stimulus.TimeAxis*Stimulus.Frequency,Stimulus.Signal);
              
               plot(current,Stimulus.Signal);
               figure;
               plot(current)
            xlabel('Period Number');
%             ylabel(axesHandles(1),'Voltage (V)');
%             ylabel(axesHandles(2),'Current (A)');
            title('Current as func. of period number');
            curLine.Color='red';
            curLine.LineWidth = 1;
            volLine.Color='cyan';
            volLine.LineWidth = 0.1;
            grid on;set(gca,'FontSize',16);  set(gca,'XMinorTick','on','YMinorTick','on');
% 
%             subplot(2,1,2);
if plotpsd
             figure,
            hermonicsVecor = t_freq/Stimulus.Frequency;
            extCur_freq = extCur_freq(hermonicsVecor<10);
            hermonicsVecor = hermonicsVecor(hermonicsVecor<10);
            plot(hermonicsVecor,extCur_freq);
            set(gca,'XTick',[0:1:10]);
            xlabel('Harmonic')
            ylabel('psd of current(s2/hz)')
            title('Power density spectrum of the current')
            set(gca,'FontSize',16);  set(gca,'XMinorTick','on','YMinorTick','on');
            grid on;
end
end
        case 'DCandWait'
            [axesHandles, ~, ~] = plotyy(Stimulus.TimeAxis,conductance,...
                                         Stimulus.TimeAxis,Stimulus.Signal);
            title('External voltage and conductance as function of time')
            xlabel('Time (sec)');
            ylabel(axesHandles(1),'Conductance (S)');
            ylabel(axesHandles(2),'External Voltage (V)');
        case 'Ramp'
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
end