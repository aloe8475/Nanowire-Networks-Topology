function Stimulus = getStimulus(Stimulus, SimulationOptions)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generates a structure with the details of an external voltage signal
% applied to the network.
%
% ARGUMENTS: 
% Stimulus - Structure containing the details of the required stimulus. It
%           must contain a field .BiasType, whose value is a string
%           specifying the type of stimulus.:
%           - 'DC' - Constant external voltage. 

%                      .Amplitude
%           - 'AC' - Sinusoidal external voltage.
%                    Required fields:
%                      .Frequency
%                      .Amplitude
%           - 'DCandWait' - A DC signal followed by a much smaller DC
%                           signal (which is meant only for probing, rather 
%                           then for stimulating).
%                           Required fields:
%                             .T 
%                             .dt
%                             .OffTime % Time at which tthe stimulus changes to AmplitudeOff
%                             .AmplitudeOn
%                             .AmplitudeOff
%           - 'Ramp' - A ramping voltage signal.
%                      Required fields:
%                        .T
%                        .dt
%                        .AmplitudeMin
%                        .AmplitudeMax
%           - 'Custom' - An arbitrary voltage signal.
%                        Required fields:
%                          .T
%                          .dt
%                          .TimeAxis
%                          .Signal
% SimulationOptions -  a struct with additional information about the simulation
%                    Required fields:
%                      .T                           (signal duration)
%                      .dt                          (duration of time-step)
%           It is assumed that the units of all input fields are sec, Hz
%           and Volt. Thus, the output fields are in sec and Volt.
% OUTPUT:
% Stimulus - Structure with the details of the time axis and with the 
%            external voltage signal. Fields:
%              .BiasType
%              .T
%              .dt
%              .TimeAxis
%              .Signal
%              .Frequency [Hz] (for AC and Sawtooth)
%
% REQUIRES:
% none
%
% USAGE:
%{
    Options.T          = 1e+1; % (sec)
    Options.dt         = 1e-3; % (sec)
    Stimulus.BiasType  = 'DC';       
    Stimulus.Amplitude = 3;    % (Volt)

    Stimulus = getStimulus(Stimulus, Options);
%}
%
% Authors:
% Ido Marcus
% Paula Sanz-Leon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    % Signal duration and timestep:
    Stimulus.T  = SimulationOptions.T;
    Stimulus.dt = SimulationOptions.dt;
    Stimulus.TimeAxis = SimulationOptions.TimeVector;
    % = (Stimulus.dt : Stimulus.dt : Stimulus.T)';
    % The first time point we can actually records happens at dt. 
    % Time=0 correspond for the values of the intial conditions
    
    % External voltage signal:
    switch Stimulus.BiasType
        case 'DOW'
%             Stimulus.Signal = Stimulus.Amplitude*ones(size(Stimulus.TimeAxis));
    T=dlmread('C:\Users\aloe8475\Documents\PhD\GitHub\CODE\Analysis\Classification\Time Serieis\Dow_Jones_Percentage_Change1.csv');
    T = Stimulus.Amplitude*T;
    %     T(T==0)=-0.5;
%     T(T==1)=Stimulus.Amplitude;
    for i=1:length(Stimulus.TimeAxis)
        Stimulus.Signal(i)=T(i);
%        Stimulus.Signal(i) =T(ceil(i/10));
    end
    Stimulus.Signal=Stimulus.Signal';
        case 'AC'
%            Stimulus.Signal = Stimulus.Amplitude*sin(2*pi*Stimulus.Frequency*Stimulus.TimeAxis);
            Stimulus.Signal = Stimulus.Amplitude*sawtooth(2*pi*Stimulus.Frequency*(Stimulus.TimeAxis-0.75/Stimulus.Frequency) , 0.5);
        case 'DCandWait'
            % Standard boxcar:
            %Stimulus.Signal = Stimulus.AmplitudeOn*ones(size(Stimulus.TimeAxis));
            % rectangular pulse train:
            Stimulus.Signal = max(Stimulus.AmplitudeOff,Stimulus.AmplitudeOn*square(1*pi*Stimulus.TimeAxis));
            % sawtooth (for I-V plots):
            %Stimulus.Signal = Stimulus.AmplitudeOn*sawtooth(2*pi*1.0*(Stimulus.TimeAxis-0.75/10.0) , 0.5);
%            % Walsh functions:
%            Stimulus.Signal = max(Stimulus.AmplitudeOff,Stimulus.AmplitudeOn*sign(cos(8*pi*Stimulus.TimeAxis).*sin(2*pi*Stimulus.TimeAxis))); % Walsh9
%            Stimulus.Signal = max(Stimulus.AmplitudeOff,Stimulus.AmplitudeOn*sign(cos(2*pi*Stimulus.TimeAxis).*cos(4*pi*Stimulus.TimeAxis)));
%            % Walsh6
%            Stimulus.Signal =
%            max(Stimulus.AmplitudeOff,Stimulus.AmplitudeOn*sign(cos(2*pi*Stimulus.TimeAxis).*sin(2*pi*Stimulus.TimeAxis)));
%            spike:
%            Stimulus.Signal = Stimulus.AmplitudeOff*ones(size(Stimulus.TimeAxis));
%            Stimulus.Signal(Stimulus.TimeAxis >= 5) = max( Stimulus.AmplitudeOff, exp( -1.5*((Stimulus.TimeAxis(Stimulus.TimeAxis >= 5)/5)- 1) )*Stimulus.AmplitudeOn );
            %
            Stimulus.Signal(Stimulus.TimeAxis > Stimulus.OffTime) = Stimulus.AmplitudeOff;
        case 'Ramp'
            Stimulus.Signal = linspace(Stimulus.AmplitudeMin, Stimulus.AmplitudeMax, length(Stimulus.TimeAxis))';
        case 'Custom'
            Stimulus.Signal = Stimulus.Signal;
        case 'Triangle'
            t  = SimulationOptions.TimeVector - SimulationOptions.dt - SimulationOptions.T/2; 
            w  = SimulationOptions.T;
            scaling = abs(Stimulus.AmplitudeMax - Stimulus.AmplitudeMin);
            Stimulus.Signal = scaling * tripuls(t, w) + Stimulus.AmplitudeMin;
    end
end