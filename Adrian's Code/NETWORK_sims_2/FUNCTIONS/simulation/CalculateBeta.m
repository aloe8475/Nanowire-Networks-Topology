function beta=CalculateBeta(Sim)

current=Sim.Data.Isource;
time=Sim.Data.time;
[fp,sp]=psd_local(full(current),full(time));
beta=fast_beta(fp,sp',50);

end