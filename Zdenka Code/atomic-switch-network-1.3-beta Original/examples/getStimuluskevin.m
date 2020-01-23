%%
% 
% 
% 
function [glo,Stimulus] = getStimuluskevin(glo,Stimulus, SimulationOptions)
% global tau;
% global ap;global c;
% global b;
% global M;global t;global dutyratio;
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
   t0=0.1*100;
    Stimulus.T  = SimulationOptions.T;
    Stimulus.dt = SimulationOptions.dt;
    Stimulus.TimeAxis = SimulationOptions.TimeVector;
    b=glo.b;
    % = (Stimulus.dt : Stimulus.dt : Stimulus.T)';
    % The first time point we can actually records happens at dt. 
    % Time=0 correspond for the values of the intial conditions
    
    % External voltage signal:
    switch Stimulus.BiasType
        case 'DC'
%             Stimulus.Signal =ap.*Stimulus.Amplitude*mackeyglass(sum(size(Stimulus.TimeAxis))-1,0,0.2,0.1,10)-b;
% Stimulus.Signal = ap.*ones(size(Stimulus.TimeAxis));
%  Stimulus.Signal(Stimulus.TimeAxis > c)=b*ones(size(Stimulus.TimeAxis(Stimulus.TimeAxis > c)));
for aa=1:glo.dutyratio-1
    for i=1:length(Stimulus.TimeAxis)
          Stimulus.Signal(aa,i) = max(0,glo.ap*(glo.M{b}(aa,max(1,ceil(i/t0)))));
    end
end
Stimulus.Signal(end+1,:)=0;
  for c=2:length(Stimulus.TimeAxis)
      if(~sum(Stimulus.Signal(:,c-1))==0&&sum(Stimulus.Signal(:,c))==0)
%           Stimulus.Signal(:,c)=1;
          
           break;
      end
  end
  glo.c=c;
  for t=2:length(Stimulus.TimeAxis)
      if(sum(Stimulus.Signal(:,t-1))==0&&~sum(Stimulus.Signal(:,t))==0)
           break;
      end
  end
  glo.t=t;
        case 'AC'
            cc=glo.ap.*ones(size(Stimulus.TimeAxis));
%             aa=1*(glo.ap.*mackeyglass(length(Stimulus.TimeAxis)+2000-1,glo.tau))-glo.th;
%             aa(1:2000)=[];
            bb=(glo.ap*sin(2*glo.a*pi*Stimulus.Frequency*Stimulus.TimeAxis));
            Stimulus.Signal=zeros(1,length(Stimulus.TimeAxis));
            Stimulus.Signal(1,:)=bb;
%             for a=1:length(SimulationOptions.ContactNodes)-1
%             Stimulus.Signal(a,:)=bb;
%             end

%  Stimulus.Signal =(glo.ap*sin(2*glo.a*pi*Stimulus.Frequency*Stimulus.TimeAxis));
%  Stimulus.Signal(end+1,:)=0;
%   Stimulus.Signal2 =(ap*sin(2*a*pi*Stimulus.Frequency*Stimulus.TimeAxis));
%             Stimulus.Signal(Stimulus.TimeAxis > 20) =(ap**sin(2*a*pi*Stimulus.Frequency*Stimulus.TimeAxis(Stimulus.TimeAxis > 20))+b);
           Stimulus.Signal(end+1,:)=0;
        case 'DCandWait'
%             Stimulus.Signal = a*ones(size(Stimulus.TimeAxis));
 %             1. sin series
%             Stimulus.Signal =max(Stimulus.AmplitudeOff,(ap*sin(2*a*pi*Stimulus.TimeAxis)+ap));
%             Stimulus.Signal(Stimulus.TimeAxis > 30) =max(Stimulus.AmplitudeOff,(ap*sin(2*a*pi*Stimulus.TimeAxis(Stimulus.TimeAxis > 30))+ap)...
%                 .*(square(b*Stimulus.TimeAxis(Stimulus.TimeAxis > 30),dutyratio)+1)/2);
%    Stimulus.Signal =max(Stimulus.AmplitudeOff,(ap*sin(2*a*pi*Stimulus.TimeAxis)+ap));
%       Stimulus.Signal =(ap*sin(2*a*pi*Stimulus.TimeAxis)+ap);
%               Stimulus.Signal=1*mackeyglass(sum(size(Stimulus.TimeAxis))-1,0,0.2,0.1,10)-0.6;
%                Stimulus.Signal(rem(Stimulus.TimeAxis,0.07)==0)=Stimulus.AmplitudeOff;
%   Stimulus.Signal(Stimulus.TimeAxis > c) =b*ones(size(Stimulus.TimeAxis(Stimulus.TimeAxis > c)));
%    Stimulus.Signal(Stimulus.TimeAxis > c+10) =max(Stimulus.AmplitudeOff,(ap*ones(size(Stimulus.TimeAxis(Stimulus.TimeAxis > c+10)))));
%    Stimulus.Signal(Stimulus.TimeAxis > c+10) =max(Stimulus.AmplitudeOff,(ap*sin(2*a*pi*Stimulus.TimeAxis(Stimulus.TimeAxis > c+10))+ap));
%    Stimulus.Signal(Stimulus.TimeAxis > c+90) =max(Stimulus.AmplitudeOff,(0.15*sin(2*a*pi*Stimulus.TimeAxis(Stimulus.TimeAxis > c+90))+0.15));
%             Stimulus.Signal(Stimulus.TimeAxis > 150) =0.006*ones(size(Stimulus.TimeAxis(Stimulus.TimeAxis > 150)));
%             Stimulus.Signal(Stimulus.TimeAxis > 200) =0.0055*ones(size(Stimulus.TimeAxis(Stimulus.TimeAxis > 200)));
%             Stimulus.Signal(Stimulus.TimeAxis > 250) =0.005*ones(size(Stimulus.TimeAxis(Stimulus.TimeAxis > 250)));
%             Stimulus.Signal(Stimulus.TimeAxis > 300) =0.0045*ones(size(Stimulus.TimeAxis(Stimulus.TimeAxis > 300)));
%             Stimulus.Signal(Stimulus.TimeAxis > 350) =0.004*ones(size(Stimulus.TimeAxis(Stimulus.TimeAxis > 350)));
%            Stimulus.Signal(Stimulus.TimeAxis > 450) =0.0035*ones(size(Stimulus.TimeAxis(Stimulus.TimeAxis > 450)));
% Stimulus.Signal=max((ap*sin(2*a*pi*Stimulus.TimeAxis)+ap).*((square(b*Stimulus.TimeAxis,dutyratio)+1)/2),0.04*sin(2*a*pi* Stimulus.TimeAxis)+0.04);
% Stimulus.Signal(Stimulus.TimeAxis>300) =2*ones(size(Stimulus.TimeAxis(Stimulus.TimeAxis>300)));
% Stimulus.Signal(Stimulus.TimeAxis>400)=max((ap*sin(2*a*pi*Stimulus.TimeAxis(Stimulus.TimeAxis>400))+ap).*((square(b*Stimulus.TimeAxis(Stimulus.TimeAxis>400),dutyratio)+1)/2),0.04*sin(2*a*pi* Stimulus.TimeAxis(Stimulus.TimeAxis>400))+0.04);
%             Stimulus.Signal(Stimulus.TimeAxis > 20)=0.01.*mackeyglass(sum(size(Stimulus.TimeAxis(Stimulus.TimeAxis>20)))-1,0,0.2,0.1,10);
% Stimulus.Signal(Stimulus.TimeAxis>10)=max((0.17*sin(2*a*pi*Stimulus.TimeAxis(Stimulus.TimeAxis>10))+0.17).*((square(b*Stimulus.TimeAxis(Stimulus.TimeAxis>10),dutyratio)+1)/2),0.04*sin(2*a*pi* Stimulus.TimeAxis(Stimulus.TimeAxis>10))+0.04);
% Stimulus.Signal(Stimulus.TimeAxis > 60)=max(Stimulus.AmplitudeOff,0.1*sin(2*a*pi* Stimulus.TimeAxis(Stimulus.TimeAxis>60))+0.1);
% Stimulus.Signal(Stimulus.TimeAxis > 100)=max(Stimulus.AmplitudeOff,0.2*sin(2*a*pi* Stimulus.TimeAxis(Stimulus.TimeAxis>100))+0.2);
             %             Stimulus.Signal(Stimulus.TimeAxis > 1000)=(ap-b-0.04)*sin(2*a*pi* Stimulus.TimeAxis(Stimulus.TimeAxis>1000))+ap-b-0.03;
%             Stimulus.Signal(Stimulus.TimeAxis > 1500)=(ap-b-0.00)*sin(2*a*pi* Stimulus.TimeAxis(Stimulus.TimeAxis>1500))+ap-b-0.00;
            % 2.pulsed v
Stimulus.Signal = glo.ap.*ones(size(Stimulus.TimeAxis));
%  Stimulus.Signal(Stimulus.TimeAxis > glo.c)=b*ones(size(Stimulus.TimeAxis(Stimulus.TimeAxis > c)));

%                          Stimulus.Signal = ap.*max(Stimulus.AmplitudeOff,Stimulus.AmplitudeOn.*square(2*a*pi*Stimulus.TimeAxis,dutyratio)...
%                  .*((square(c.*Stimulus.TimeAxis)+1)./2));
% 3.trangle wave
%              Stimulus.Signal = Stimulus.AmplitudeOn*sawtooth(a*pi*1.0*(Stimulus.TimeAxis-0.75/10.0) , 0.5);
% 4.dc
%             Stimulus.Signal = 0.02*ones(size(Stimulus.TimeAxis));
%5.walsh
% Stimulus.Signal = max(Stimulus.AmplitudeOff,0*sign(cos(a*pi*Stimulus.TimeAxis).*sin(2*pi*Stimulus.TimeAxis)));
%       spike:
%           Stimulus.Signal = Stimulus.AmplitudeOff*ones(size(Stimulus.TimeAxis));
%           Stimulus.Signal(Stimulus.TimeAxis >= 5) = max( Stimulus.AmplitudeOff, 0.6*exp( -1.5*((Stimulus.TimeAxis(Stimulus.TimeAxis >= 5)/5)- 1) )*Stimulus.AmplitudeOn );
%             
        
        case 'Ramp'
            
%             Stimulus.Signal = linspace(Stimulus.AmplitudeMin, Stimulus.AmplitudeMax, length(Stimulus.TimeAxis))';
        case 'Custom'
            Stimulus.Signal = Stimulus.Signal;
        case 'Triangle'
            t  = SimulationOptions.TimeVector - SimulationOptions.dt - SimulationOptions.T/2; 
            w  = SimulationOptions.T;
            scaling = abs(Stimulus.AmplitudeMax - Stimulus.AmplitudeMin);
            Stimulus.Signal = scaling * tripuls(t, w) + Stimulus.AmplitudeMin;
    end
end