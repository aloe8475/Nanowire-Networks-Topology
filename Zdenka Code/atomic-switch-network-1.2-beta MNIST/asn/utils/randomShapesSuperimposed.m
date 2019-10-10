function outputData = randomShapesSuperimposed(iter,shape,withPlot)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function superimposes squares\Gaussian of random position, width and 
% height, computes the PSD of the resulting random signal and returns the 
% slope of the PSD as it appears in a log-log plot.
%
% ARGUMENTS: 
% iter - number of iterations (each iteration generates and analysis a new
%        random signal).
% shape - 'square' \ 'Gaussian' (no default, must be specified).
% withPlot - true if signal and PSD are to be plotted.
%
% OUTPUT:
% outputData - a structutre containing the result. fields:
%                .slope - a 1Xiter vector with the above mentioned slopes.
%                .signal - a 1Xiter cell array with the generated random 
%                          signals.
%                .PSD - a 1Xiter cell array with the PSD of the random 
%                       signals.
%                .timeline - x axis for the signals.
%                .freqline - x axis for the PSDs.
%
% REQUIRES:
% fourier_analysis
%
% Authors:
% Ido Marcus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialize
% Allocate space for results:
outputData.slope  = zeros(1,iter);
outputData.signal = cell(1,iter);
outputData.PSD    = cell(1,iter);

% Signal parameters:
dt = 1e-4;        % Time-step          (sec)
T  = 1e+3;        % Duration of signal (sec)

% Time axis:
t  = 0:dt:T;      % time-axis          (sec)    
if mod(length(t),2) == 1    % a real signal of even length will have a 'Hermitian' Fourier transform
    t = t(1:end-1);
end

% Shape parameters
n  = 1e1;                               % Mean number of shapes per second of signal (1/sec, Hz)
%durationMean = ones(iter,1)*1e+0;      % Mean duration of a single shape            (sec)
durationMean = [50,10,2,0.4];           % Mean duration of a single shape            (sec)
durationStd  = ones(iter,1)*1e-1;       % Standard deviation of the mean of duration (sec)

% Analysis parameters
noiseThreshold = 1e1; % Lowest frequency to be used in the linear fitting of the log-log PSD (Hz)

%% Run:
close all;
figure;
N = T*n;
for j=1:iter
    % Progress bar:
    fprintf('\nIteration %d\\%d ', j, iter);
    
    % Shapes random parameters:
    location = randi(length(t),N,1);                                        % When shapes start                (index)
    height   = rand(N,1);                                                   % Height of shapes                 (arbitrary units)
    width = 2*floor(normrnd(durationMean(j)/dt,durationStd(j)/dt,N,1)/2)+1; % Duration of shapes (odd number)  (index)  
    
    % Superimpose squares:
    signal = zeros(1,2*length(t)); % Signal is too long on purpose, as it allows periodic boundary conditions.
    for i=1:N
        progressBar(i,N);
        
        currIndices = location(i):location(i)+width(i)-1;
        switch shape
            case 'square'
                currShape = height(i)*ones(1,width(i));
            case 'Gaussian'
                currShape = Gaussian(-floor(width(i)/2):floor(width(i)/2),width(i)/8,height(i));
        end
        signal(currIndices) = signal(currIndices) + currShape;
    end
    signal = signal(1:end/2)+signal(end/2+1:end); % Apply periodic boundary conditions.
    %signal = signal-mean(signal);                 % Get rid of zero-frequency component.

    % Fourier trnasform the result:
    [freqVector, freqSignal] = fourier_analysis(t, signal);

    % Linear fit of PSD in log-log scale:
    coefficients = polyfit(log10(freqVector(freqVector>noiseThreshold)),log10(freqSignal(freqVector>noiseThreshold)),1);
    fitPSD = 10^coefficients(2)*freqVector.^(coefficients(1));

    % Save results:
    if j == 1
        outputData.timeline = t;
        outputData.freqline = freqVector;
    end
    outputData.signal{j} = signal;
    outputData.PSD{j}    = freqSignal;
    outputData.slope(j)  = -coefficients(1);
    
    % Plot results if asked to:
    if withPlot      
        subplot(iter,2,(j-1)*2+1)
        plot(t,signal);
        switch shape
            case 'square'
                title(sprintf('Superimposed squares of average length = %d sec', durationMean(j)));
            case 'Gaussian'
                title(sprintf('Superimposed Gaussians of average length = %d sec', durationMean(j)));
        end
        xlabel('Time (sec)');
        ylabel('Height (au)');
        grid on;

        subplot(iter,2,(j-1)*2+2)
        loglog(freqVector(freqVector<noiseThreshold),freqSignal(freqVector<noiseThreshold),'Color','b');
        hold on;
        loglog(freqVector(freqVector>=noiseThreshold),freqSignal(freqVector>=noiseThreshold),'Color','g');
        hold on;
        loglog(freqVector,fitPSD);
        title('Power density spectrum of signal');
        xlabel('Frequency (Hz)');
        ylabel('Power density');
        legend({'PSD','PSD used for fitting', 'Linear fit'});
        xlim([1/T,1/dt/2]);
        ylim([1e-11,1e1]);
        set(gca,'Ytick',10.^(-12:2:12));
        text(1e0,1e-2,sprintf('\\beta=%.2f', -coefficients(1)),'Color','r','FontSize',18);
        grid on;
        
%         if strcmp(shape,'square')
%             [interval_hist_bin, interval_hist_values] = find_stability_histogram(signal',dt);
%             
%             subplot(2,2,3);
%             plot(interval_hist_bin, interval_hist_values,'*');
%             title('Histogram of lengths of intervals between changes in signal');
%             xlabel('Interval length (sec)')
%             ylabel('Probability density function');
%             grid on;
% 
%             subplot(2,2,4);
%             semilogy(interval_hist_bin, interval_hist_values,'*');
%             title('Histogram of lengths of intervals between changes in signal (semilog)');
%             xlabel('Interval length (sec)')
%             ylabel('Probability density function');
%             grid on;
%         end
    end
end