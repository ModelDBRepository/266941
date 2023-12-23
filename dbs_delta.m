function I_dbs=dbs_delta(fdbs,dbs_duration,dev,nSim,Fs,cut)

%This is to define Dirac delta pulses, no membrane current but straight dirac delta pulses that reach PNs:
Tdbs=Fs/fdbs;
dbs=1:round(Tdbs):dbs_duration;
I_dbs_full=zeros(1,dbs_duration);
I_dbs_full(dbs)=1; %Amplitude of DBS

if dev==1
    I_dbs=I_dbs_full;
else
    I_dbs=[zeros(1,cut),zeros(1,(nSim-cut)/dev),I_dbs_full,zeros(1,(nSim-cut)/dev)]; %extracellular dbs pulsing
end
