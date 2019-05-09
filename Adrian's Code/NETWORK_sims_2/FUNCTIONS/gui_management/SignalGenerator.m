function [t,v,open]=SignalGenerator(settings)

Ttime=settings.Time;
Tstep=settings.Step;
Vmax=settings.Vmax;
Vmin=settings.Vmin;
NoC=settings.NoC;
Frequency=settings.Frequency;
if ~isequal(Frequency,0)
    NoC=round(Ttime*Frequency);
end
duty=settings.Duty/100;
type=settings.SigType;
TStartOff=settings.TStartOff;
SetFreq=settings.SetFreq;
VSet=settings.VSet;
ReadFreq=settings.ReadFreq;
VRead=settings.VRead;
VStartOff=settings.VStartOff;
NoP=Ttime/Tstep;
t=0:Tstep:Ttime-Tstep;
switch type
    case 'Constant'
        v=ones(1,NoP)*Vmax;
    case 'Triangular(IV)'
        T=max(t)/NoC;
        v=zeros(1,NoP);
        for i=1:length(t)
            phase=rem(t(i),T);
            v(i)=triangulo(phase,Vmax,Vmin,T);
        end
        
    case 'Square'
        T=max(t)/NoC;
        v=zeros(1,NoP);
        for i=1:length(t)
            phase=rem(t(i),T);
            v(i)=cuadrado(phase,Vmax,Vmin,T,duty);
        end
        v(end)=Vmin;
    case 'SquarePWM'
        T=max(t)/NoC;
        v=zeros(1,NoP);
        for i=1:length(t)
            phase=rem(t(i),T);
            v(i)=cuadrado(phase,1,0,T,duty);
        end
        v(end)=0;
        t1=t(v==1);T1=1./SetFreq;
        t2=t(v==0);T2=1./ReadFreq;
        SetDuty=settings.SetDuty/100;ReadDuty=settings.ReadDuty/100;
        v1=zeros(1,length(t1));
        v2=zeros(1,length(t2));
        for i=1:length(t1)
            phase=rem(t1(i),T1);
            v1(i)=cuadrado(phase,VSet,Vmin,T1,SetDuty);
        end
        for i=1:length(t2)
            phase=rem(t2(i),T2);
            v2(i)=cuadrado(phase,VRead,Vmin,T2,ReadDuty);
        end
        t=[t1 t2];
        v=[v1 v2];
        [t, idx]=sort(t);
        v=v(idx);
        
        
        
        
        
        
end
if TStartOff>0
    ta=0:Tstep:TStartOff-Tstep;
    t=t+TStartOff;
    va=ones(1,length(ta)).*VStartOff;
    t=[ta t];
    v=[va v];
end
t=t(1:NoP);
v=v(1:NoP); 
open=zeros(1,length(t));
if isequal(settings.OpenCheck,1)
    open(v==0)=1;
end



    function vval=triangulo(tt,Vmax,Vmin,T)
        if tt<=T/2
            vval=(Vmax)*2*tt/T+Vmin;
        else
            vval=-(Vmax)*(2*tt/T-2)+Vmin;
        end
    end
    function vval=cuadrado(tt,Vmax,Vmin,T,duty)
        if tt<=T*duty
            vval=Vmax;
        else
            vval=Vmin;
        end
    end
end