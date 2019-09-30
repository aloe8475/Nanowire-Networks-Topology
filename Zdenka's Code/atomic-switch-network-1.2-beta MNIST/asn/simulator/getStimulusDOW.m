
function Stimulus = getStimulus(a,Stimulus, SimulationOptions)
global dutyratio;
global ap;global c;
global b;
global M;global t;
        
    % Signal duration and timestep:
   t0=0.1*100;
    Stimulus.T  = SimulationOptions.T;
    Stimulus.dt = SimulationOptions.dt;
    Stimulus.TimeAxis = SimulationOptions.TimeVector;
    % = (Stimulus.dt : Stimulus.dt : Stimulus.T)';
    % The first time point we can actually records happens at dt. 
    % Time=0 correspond for the values of the intial conditions
    
    % External voltage signal:
    switch Stimulus.BiasType
        case 'DC'
%             Stimulus.Signal =ap.*Stimulus.Amplitude*mackeyglass(sum(size(Stimulus.TimeAxis))-1,0,0.2,0.1,10)-b;
% Stimulus.Signal = ap.*ones(size(Stimulus.TimeAxis));
%  Stimulus.Signal(Stimulus.TimeAxis > c)=b*ones(size(Stimulus.TimeAxis(Stimulus.TimeAxis > c)));
T=dlmread('C:\Users\aloe8475\Documents\PhD\GitHub\CODE\Analysis\Classification\Time Serieis\Dow_Jones_Percentage_Change.csv');
    for i=1:length(Stimulus.TimeAxis)
          Stimulus.Signal(i) =T(i);
    end
        case 'AC'

            Stimulus.Signal=2*(ap.*mackeyglass(length(Stimulus.TimeAxis)-1+1000,17))+1;Stimulus.Signal=Stimulus.Signal(1001:end);
%  Stimulus.Signal =(ap*sin(2*a*pi*Stimulus.Frequency*Stimulus.TimeAxis));
%   Stimulus.Signal2 =(ap*sin(2*a*pi*Stimulus.Frequency*Stimulus.TimeAxis));
%             Stimulus.Signal(Stimulus.TimeAxis > 20) =(ap**sin(2*a*pi*Stimulus.Frequency*Stimulus.TimeAxis(Stimulus.TimeAxis > 20))+b);
           
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
                        Stimulus.Signal = ap*max(Stimulus.AmplitudeOff,Stimulus.AmplitudeOn*square(2*a*pi*Stimulus.TimeAxis,dutyratio));
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