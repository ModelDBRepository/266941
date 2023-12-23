clearvars
dd=datestr(datetime); dd=[dd(1:14),dd(16:17)]
mkdir(dd);

%Number of trials:
Nt=1;
% TURN DBS on OR off:
% fidD=5*67; %DBS ON with the synaptic fidelity fidD, for dbs carriers (to be used to invade layer D)
fidD=0;  %DBS OFF
%Set the number of affected PNs in D by DBS
nh=0.1;     %percentage of PNs that are hyperdirect

%%% SIMULATION TIME AND STEP, Time delays~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sim_time=10;          %Simulation time in seconds (must be a multiplacative of 3 under PD+DBS condition)
T=(sim_time+1)*1000; %Simulation time in ms with 1 extra second to reach the steady state and trash later
dt=.1;               %time span and step (of ms)
Fs=1000/dt;          %sampling frequency in Hz
nSim=round(T/dt);    %number of simulation steps
chop_till=1*Fs;      %Cut the first 1 seconds of the simulation
td_L=8;      %8      %time delay between the layers in corticex and nuclei in thalamus (ms)
td_wL=1;     %1      %time delay within a structure (ms)
td_TC=15;    %25     %time delay between thalamus and cortex (ms) (transmission time delay)
td_CT=20;            %time delay between cortex and thalamus (ms) (transmission time delay)  
td_syn=1;    %x      %Synaptic transmission delay (fixed for all synapses in the TCM)
if td_TC>=td_CT
tVec=td_TC+td_syn+1:nSim;     %time vector
end
if td_TC<td_CT
tVec=td_CT+td_syn+1:nSim;     %time vector
end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ECoG=zeros(Nt,nSim-chop_till);
LFP=zeros(Nt,nSim-chop_till);
X=zeros(6,nSim-chop_till,Nt); 

for ij=1:Nt
    
display(ij,'Trial')

%Number of Excitatory  & Inhiibitory Neurons in the structures, nEs is # of neurons in S, nEm is # of neurons in M, ...
% nEs=10; nEm=10; nEd=10; nErel=10; nINs=10; nIret=4;
% nEs=50; nEm=50; nEd=50; nErel=50; nINs=50; nIret=20;
nEs=100; nEm=100; nEd=100; nErel=100; nINs=100; nIret=40;
% nEs=200; nEm=200; nEd=200; nErel=200; nINs=200; nIret=80;
% nEs=500; nEm=500; nEd=500; nErel=500; nINs=500; nIret=200;
% nEs=1000; nEm=1000; nEd=1000; nErel=1000; nINs=1000; nIret=400;
% nEs=2000; nEm=2000; nEd=2000; nErel=2000; nINs=2000; nIret=800;
% nEs=600; nEm=600; nEd=600; nErel=600; nINs=600; nIret=240;
n_tot=nEs+nEm+nEm+nINs+nErel+nIret;

%Impact of DBS on the other cortical structures via D PNs axons:
fidCI=abs(1*fidD);   %the synaptic fidelity, for dbs carriers (to be used to invade CIs)
fidM=abs(0*fidD);    %the synaptic fidelity, for dbs carriers (to be used to invade layer M)
fidS=abs(1*fidD);    %the synaptic fidelity, for dbs carriers (to be used to invade layer S)
fidR=abs(1*fidD);    %the synaptic fidelity, for dbs carriers (to be used to invade layer TCR)
fidN=abs(1*fidD);    %the synaptic fidelity, for dbs carriers (to be used to invade layer TRN)
nCI=1; nS=1; nR=1; nN=1;
n_conn_CI=nCI*nh*nINs;  % percentage of CI neurons that have synaptic contact with hyperdirect neurons axon arbors
n_conn_S=nS*nh*nEs;    % percentage of S neurons that have synaptic contact with hyperdirect neurons axon arbors
n_conn_M=0*nh*nEm;    % percentage of M neurons that have synaptic contact with hyperdirect neurons axon arbors
n_conn_R=nR*nh*nErel;  % percentage of R neurons that have synaptic contact with hyperdirect neurons axon arbors
n_conn_N=nN*nh*nIret;  % percentage of N neurons that have synaptic contact with hyperdirect neurons axon arbors

n_hyp=nEd*nh; %number of hyperdirect neurons          

%Distribution of neurons in each structure                      
nE1s=.5*nEs; nE2s=.5*nEs;    nINs1=.5*nINs; nINs2=.5*nINs;
nE1m=1*nEm; nE2m= 0*nEm;     nIret1=.5*nIret; nIret2=.5*nIret;
nE1d=0.7*nEd; nE2d=0.3*nEd;   nErel1=.7*nErel; nErel2=.3*nErel;
%Types of Excitatory and Inhibitory neurons in each structure (according to
%the vectors below thw NEURONS PARAMETERS section
 E1s=1;           E2s=2;     I1s=4;             I2s=5; 
 E1m=1;           E2m=1;     I1ret=8;           I2ret=8;
 E1d=1;           E2d=2;     E1rel=6;           E2rel=6;
% NEURONS PARAMETERS
% 1)Rs    2) IB   3)CH    4)FS   5)LTS   6)Rel(TC)  7)Rel(CH)  8)Ret
a=[0.02,  0.02,   0.05,   0.1,   0.02,   0.02,      0.02,      0.02]; 
b=[0.2,   0.2,    0.2,    0.2,   0.25,   0.25,      0.25,      0.25]; 
c=[-65,   -55,    -50,    -65,   -65,    -65,       -65,        -65];
d=[8,     4,      2,      2,      2,     0.05,      0.05,      2.05]; 

%%%connectivity factor
% fac_PD=0.25; %for 10 neurons
fac_N=2.5; fac_PD=5; %for 100 neurons 
% fac_PD=.3; %for 500 neurons
% fac_PD=.1;  %for 1000 neurons

%%%Choose the .m file that contains the coupling matrix:
% PDone    %Contains coupling strengths to mimic high-beta state (PD state)
NormalOne  %Contians coupling strengths to mimic low-beta state (non-PD state)
%%%%%%%%%%%%%%%%%%%%%%%%%% 
if ij==Nt
%%%Construct the normalized mean synaptic weights, (Fig. 1B and C in the
%%%IEEE conference paper)
%1=Layer S, 2=Layer M, 3=Layer D, 4=INs, 5=TRN, 6=TCR
Z(1,1)=mean(W_EEs); Z(2,2)=mean(W_EEm); Z(3,3)=mean(W_EEd); z(4,4)=-mean(W_IIins); Z(6,6)=mean(W_IIret); Z(5,5)=mean(W_EErel);
Z(2,1)=mean(W_EEsm); Z(3,1)=mean(W_EEsd); Z(4,1)=mean(W_EIsINs); Z(6,1)=mean(W_EIsRet); Z(5,1)=mean(W_EEsRel);
Z(1,2)=mean(W_EEms); Z(3,2)=mean(W_EEmd); Z(4,2)=mean(W_EImINs); Z(6,2)=mean(W_EImRet); Z(5,2)=mean(W_EEmRel);
Z(1,3)=mean(W_EEds); Z(2,3)=mean(W_EEdm); Z(4,3)=mean(W_EIdINs); Z(6,3)=mean(W_EIdRet); Z(5,3)=mean(W_EEdRel);
Z(1,4)=mean(W_IE_INs_s); Z(2,4)=mean(W_IE_INs_m); Z(3,4)=mean(W_IE_INs_d); Z(6,4)=mean(W_II_INs_Ret); Z(5,4)=mean(W_IE_INs_Rel);
Z(1,6)=mean(W_IE_Ret_s); Z(2,6)=mean(W_IE_Ret_m); Z(3,6)=mean(W_IE_Ret_d); Z(4,6)=mean(W_II_Ret_INs); Z(5,6)=mean(W_IE_Ret_Rel);
Z(1,5)=mean(W_EERels); Z(2,5)=mean(W_EERelm); Z(3,5)=mean(W_EEReld); Z(4,5)=mean(W_EIRelINs); Z(6,5)=mean(W_EIRelRet);
ZZ=Z; %/max(max(Z));
Zx=Z/max(max(Z));
figure; imagesc(ZZ); grid on; cmap=colormap; 
ax=gca; ax.XTick=1.5:6.5; ax.YTick=1.5:6.5;  ax.CLim=[-1 1]; col=colorbar; col.Limits=[-1 1];  ax.GridColor='k'; ax.LineWidth=1;
figdest=fullfile(dd,'TCM_Coupling_Matrix_PD');
savefig(figdest); saveas(gca,figdest,'jpeg'); saveas(gca,figdest,'svg')
close
figure; imagesc(Zx); grid on; cmap=colormap; 
ax=gca; ax.XTick=1.5:6.5; ax.YTick=1.5:6.5;  ax.CLim=[-1 1]; col=colorbar; col.Limits=[-1 1];  ax.GridColor='k'; ax.LineWidth=1;
figdest=fullfile(dd,'TCM_Coupling_Matrix_PD_Normalized');
savefig(figdest); saveas(gca,figdest,'jpeg'); saveas(gca,figdest,'svg')
close
end
clear Z Zx
% %NO CELL VARIATION
% aEs=[a(E1s)*ones(nE1s,1);a(E2s)*ones(nE2s,1)]; bEs=[b(E1s)*ones(nE1s,1);b(E2s)*ones(nE2s,1)]; cEs=[c(E1s)*ones(nE1s,1);c(E2s)*ones(nE2s,1)]; dEs=[d(E1s)*ones(nE1s,1);d(E2s)*ones(nE2s,1)]; 
% aIs=[a(I1s)*ones(nINs1,1);a(I2s)*ones(nINs2,1)]; bIs=[b(I1s)*ones(nINs1,1);b(I2s)*ones(nINs2,1)]; cIs=[c(I1s)*ones(nINs1,1);c(I2s)*ones(nINs2,1)]; dIs=[d(I1s)*ones(nINs1,1);d(I2s)*ones(nINs2,1)]; 
% aEm=[a(E1m)*ones(nE1m,1);a(E2m)*ones(nE2m,1)]; bEm=[b(E1m)*ones(nE1m,1);b(E2m)*ones(nE2m,1)]; cEm=[c(E1m)*ones(nE1m,1);c(E2m)*ones(nE2m,1)]; dEm=[d(E1m)*ones(nE1m,1);d(E2m)*ones(nE2m,1)]; 
% aIret=[a(I1ret)*ones(nIret1,1);a(I2ret)*ones(nIret2,1)]; bIret=[b(I1ret)*ones(nIret1,1);b(I2ret)*ones(nIret2,1)]; cIret=[c(I1ret)*ones(nIret1,1);c(I2ret)*ones(nIret2,1)]; dIret=[d(I1ret)*ones(nIret1,1);d(I2ret)*ones(nIret2,1)]; 
% aEd=[a(E1d)*ones(nE1d,1);a(E2d)*ones(nE2d,1)]; bEd=[b(E1d)*ones(nE1d,1);b(E2d)*ones(nE2d,1)]; cEd=[c(E1d)*ones(nE1d,1);c(E2d)*ones(nE2d,1)]; dEd=[d(E1d)*ones(nE1d,1);d(E2d)*ones(nE2d,1)]; 
% aErel=[a(E1rel)*ones(nErel1,1);a(E2rel)*ones(nErel2,1)]; bErel=[b(E1rel)*ones(nErel1,1);b(E2rel)*ones(nErel2,1)]; cErel=[c(E1rel)*ones(nErel1,1);c(E2rel)*ones(nErel2,1)]; dErel=[d(E1rel)*ones(nErel1,1);d(E2rel)*ones(nErel2,1)]; 

%CELL VARIATIONs
% Make all cells slightliy different from each other based on the algorithm in Izhikevich 2003 IEEE paper 
re1s=rand(nE1s,1); re2s=rand(nE2s,1);  ri1s=rand(nINs1,1); ri2s=rand(nINs2,1);     
aEs=[a(E1s)*ones(nE1s,1);a(E2s)*ones(nE2s,1)];  aIs=[a(I1s)+0.008*ri1s;a(I2s)+0.008*ri2s];
bEs=[b(E1s)*ones(nE1s,1);b(E2s)*ones(nE2s,1)];  bIs=[b(I1s)-0.005*ri1s;b(I2s)-0.005*ri2s];
cEs=[c(E1s)+15*re1s.^2;c(E2s)+15*re2s.^2];      cIs=[c(I1s)*ones(nINs1,1);c(I2s)*ones(nINs2,1)]; 
dEs=[d(E1s)-.6*re1s.^2;d(E2s)-0.6*re2s.^2];     dIs=[d(I1s)*ones(nINs1,1);d(I2s)*ones(nINs2,1)]; 
%%%
re1m=rand(nE1m,1); re2m=rand(nE2m,1);  ri1m=rand(nIret1,1); ri2m=rand(nIret2,1);     
aEm=[a(E1m)*ones(nE1m,1);a(E2m)*ones(nE2m,1)];  aIret=[a(I1ret)+0.008*ri1m;a(I2ret)+0.008*ri2m];
bEm=[b(E1m)*ones(nE1m,1);b(E2m)*ones(nE2m,1)];  bIret=[b(I1ret)-0.005*ri1m;b(I2ret)-0.005*ri2m];
cEm=[c(E1m)+15*re1m.^2;c(E2m)+15*re2m.^2];      cIret=[c(I1ret)*ones(nIret1,1);c(I2ret)*ones(nIret2,1)]; 
dEm=[d(E1m)-.6*re1m.^2;d(E2m)-0.6*re2m.^2];     dIret=[d(I1ret)*ones(nIret1,1);d(I2ret)*ones(nIret2,1)]; 
%%%
re1d=rand(nE1d,1); re2d=rand(nE2d,1);  re1rel=rand(nErel1,1); re2rel=rand(nErel2,1);     
aEd=[a(E1d)*ones(nE1d,1);a(E2d)*ones(nE2d,1)];  aErel=[a(E1rel)+0.008*re1rel;a(E2rel)+0.008*re2rel];
bEd=[b(E1d)*ones(nE1d,1);b(E2d)*ones(nE2d,1)];  bErel=[b(E1rel)-0.005*re1rel;b(E2rel)-0.005*re2rel];
cEd=[c(E1d)+15*re1d.^2;c(E2d)+15*re2d.^2];      cErel=[c(E1rel)*ones(nErel1,1);c(E2rel)*ones(nErel2,1)]; 
dEd=[d(E1d)-.6*re1d.^2;d(E2d)-0.6*re2d.^2];     dErel=[d(E1rel)*ones(nErel1,1);d(E2rel)*ones(nErel2,1)]; 

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Noise terms
s2=1.5;  cn=1;  %additive white Gaussian noise strength   
e2=0.5;         %threshold white Gaussian noise strength
z2=.00;         %additive pink noise strength
kisiSE=s2*[randn(nEs,1*Fs),cn*randn(nEs,nSim-1*Fs)];       kisiSI=s2*[randn(nINs,1*Fs),cn*randn(nINs,nSim-1*Fs)];
kisiME=s2*[randn(nEm,1*Fs),cn*randn(nEm,nSim-1*Fs)];       kisiIret=s2*[randn(nIret,1*Fs),cn*randn(nIret,nSim-1*Fs)];
kisiDE=s2*[randn(nEd,1*Fs),cn*randn(nEd,nSim-1*Fs)];       kisiErel=s2*[randn(nErel,1*Fs),cn*randn(nErel,nSim-1*Fs)];
zetaSE=e2*[randn(nEs,1*Fs),cn*randn(nEs,nSim-1*Fs)];       zetaSI=e2*[randn(nINs,1*Fs),cn*randn(nINs,nSim-1*Fs)];
zetaME=e2*[randn(nEm,1*Fs),cn*randn(nEm,nSim-1*Fs)];       zetaIret=e2*[randn(nIret,1*Fs),cn*randn(nIret,nSim-1*Fs)];
zetaDE=e2*[randn(nEd,1*Fs),cn*randn(nEd,nSim-1*Fs)];       zetaErel=e2*[randn(nErel,1*Fs),cn*randn(nErel,nSim-1*Fs)];
pnSE=z2*[pinknoise(nEs,1*Fs),cn*pinknoise(nEs,nSim-1*Fs)]; pnSI=z2*[pinknoise(nINs,1*Fs),cn*pinknoise(nINs,nSim-1*Fs)];  
pnME=z2*[pinknoise(nEm,1*Fs),cn*pinknoise(nEm,nSim-1*Fs)]; pnIret=z2*[pinknoise(nIret,1*Fs),cn*pinknoise(nIret,nSim-1*Fs)];
pnDE=z2*[pinknoise(nEd,1*Fs),cn*pinknoise(nEd,nSim-1*Fs)]; pnErel=z2*[pinknoise(nErel,1*Fs),cn*pinknoise(nErel,nSim-1*Fs)];
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Idc_tune=+0.1;
% Bias currents (Subthreshold CTX and Suprethreshold THM)
Idc=[3.5,3.6,3.5,3.8,0.4,0.6,0.5,0.6]+Idc_tune*ones(1,8);

IdcS_E=[Idc(E1s)*ones(nE1s,1);Idc(E2s)*ones(nE2s,1)]; Idc_INs=[Idc(I1s)*ones(nINs1,1);Idc(I2s)*ones(nINs1,1)];
IdcM_E=[Idc(E1m)*ones(nE1m,1);Idc(E2m)*ones(nE2m,1)]; Idc_Ret=[Idc(I1ret)*ones(nIret1,1);Idc(I2ret)*ones(nIret2,1)];
IdcD_E=[Idc(E1d)*ones(nE1d,1);Idc(E2d)*ones(nE2d,1)]; Idc_Rel=[Idc(E1rel)*ones(nErel1,1);Idc(E2rel)*ones(nErel2,1)];
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%POISSONIAN background activity (assuming they invade primary motor coretex
%from premotor and supplimentary motor areas)
%Poissonian postsynaptic input to the E and I neurons for all layers
w_ps=1.0;
if w_ps==0
    I_ps=zeros(6,2,nSim);
else
    I_ps=zeros(6,2,nSim);
W_ps=w_ps*rand(6,2); 
for L=1:6
fr=20+2*randn; % Poissonian firing frequency from other parts of the brain
[spikess,tsp]=poissonSpikeGen(fr,T/1000,1,dt/1000);
tps=find(spikess==1);   
I_ps(L,1,:)=W_ps(L,1)*TMsynE(tps,nSim,td_syn,dt);
I_ps(L,2,:)=W_ps(L,2)*TMsynI(tps,nSim,td_syn,dt);
end
end
%%%%%%%DBS%%%%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
I_dbs=zeros(2,nSim); 
dev=1; %devide the total simulation time in dev sections

if fidD~=0
    dev=3;
fdbs=130;
if dev==1 %for DBS on all the time
dbs_duration=nSim;
Tdbs=Fs/fdbs;
dbs=1:round(Tdbs):dbs_duration;
I_dbs_full=zeros(1,dbs_duration);
I_dbs_full(dbs)=.02; %Amplitude of DBS
I_dbs_pre=I_dbs_full;
else
dbs_duration=(nSim-chop_till)/dev; %in seconds
I_dbs_pre=1*dbs_delta(fdbs,dbs_duration,dev,nSim,Fs,chop_till); %multiplied by 10 to make the transmembrane voltage about 80 mv.
end
%Postsynaptic DBS pulses (intra-axonal)
I_dbs_post=TMsynE_dbs(I_dbs_pre,nSim,td_syn,dt);

I_dbs(1,:)=I_dbs_pre;  
I_dbs(2,:)=I_dbs_post; 

figure; plot(I_dbs(1,:),'LineWidth',1); 
figdest=fullfile(dd,'Stimulus'); 
savefig(figdest); 
% saveas(gca,figdest,'jpeg'); saveas(gca,figdest,'svg')
close
figure; plot(I_dbs(1,:),'LineWidth',1); xlim([((nSim-10000)/dev)+10000 ((nSim-10000)/dev)+10500])
figdest=fullfile(dd,'Stimulus_zoom'); 
savefig(figdest); 
% saveas(gca,figdest,'jpeg'); 
% saveas(gca,figdest,'svg')
close
figure; plot(I_dbs(2,:),'LineWidth',1); 
figdest=fullfile(dd,'Postsynaptic_Stimulus'); 
savefig(figdest); 
% saveas(gca,figdest,'jpeg'); saveas(gca,figdest,'svg')
close
figure; plot(I_dbs(2,:),'LineWidth',1); xlim([((nSim-10000)/dev)+10000 ((nSim-10000)/dev)+10500])
figdest=fullfile(dd,'Postsynaptic_Stimulus_zoom'); 
savefig(figdest); 
% saveas(gca,figdest,'jpeg'); saveas(gca,figdest,'svg')
close
end
%%%%%%%%%%% Chnage the DBS strength ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% % for fid=2200:10:3100
% for fid=fid
%     fid
%%%%%%%%%%%%%%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Initial values for the simulation
vr=-65;
% vEs=vr*ones(1,nEs);  uEs=0*vEs;  
% vIs=vr*ones(1,nINs);  uIs=0*vIs;    
% vEm=vr*ones(1,nEm);  uEm=0*vEm;    
% vEd=vr*ones(1,nEd);  uEd=0*vEd;    
% vRet=vr*ones(1,nIret);  uRet=0*vRet;    
% vRel=vr*ones(1,nErel);  uRel=0*vRel;    
% %write v(t) and u(t) on text files
% make_v_u_textfiles
vEs=vr*ones(nEs,nSim); uEs=0*vEs; vIs=vr*ones(nEs,nSim); uIs=0*vIs;
vEm=vr*ones(nEm,nSim); uEm=0*vEm; vRet=vr*ones(nIret,nSim); uRet=0*vRet; 
vEd=vr*ones(nEd,nSim); uEd=0*vEd; vRel=vr*ones(nErel,nSim); uRel=0*vRel;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
EPSCs=zeros(1,nSim); IPSC_INs=zeros(1,nSim);
EPSCm=zeros(1,nSim); IPSC_ret=zeros(1,nSim);
EPSCd=zeros(1,nSim); EPSC_rel=zeros(1,nSim);
EPSCdF=zeros(1,nSim); EPSC_relD=zeros(1,nSim);
EPSCdD=zeros(1,nSim);

%SYNAPSE Initial values
rEs=zeros(3,1);   rEm=zeros(3,1);   rEd=zeros(3,1);  rINs=zeros(3,1);  rIret=zeros(3,1);  rErel=zeros(3,1);
xEs=ones(3,1);    xEm=ones(3,1);    xEd=ones(3,1);   xINs=ones(3,1);   xIret=ones(3,1);  xErel=ones(3,1);
IsEs=zeros(3,1);  IsEm=zeros(3,1);  IsEd=zeros(3,1); IsINs=zeros(3,1); IsIret=zeros(3,1); IsErel=zeros(3,1);
rEdF=zeros(3,1);  xEdF=ones(3,1);   IsEdF=zeros(3,1); rErelD=zeros(3,1);  xErelD=ones(3,1);   IsErelD=zeros(3,1);

%%%%% RUN THE SIMULATION %%%%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if fidD==0
    ptic('\n*** Run the simulation, numerical solution of the TCM (PD)...\n ');
else
    ptic(['\n*** Run the simulation, numerical solution of the TCM (PD+DBS@',num2str(fdbs),' Hz)...\n ']);
end
%This is coded such that we numerically solve all cells at one instance of time and then the same in the next time step and so forth
for i=tVec
%     t=i/Fs
    %Thalamic Reticular Nucleus (TRN) cells        
[vRet(:,i+1),uRet(:,i+1),rI,xI,IsI,IPSC_ret(i+1)]=...
    thm_ret(aIret,bIret,cIret,dIret,nIret,vRet(:,i),uRet(:,i),rIret,xIret,IsIret,... 
                IPSC_ret(i-td_wL-td_syn),EPSCs(i-td_CT-td_syn),EPSCm(i-td_CT-td_syn),EPSCdF(i-td_CT-td_syn),IPSC_INs(i-td_CT-td_syn),EPSC_rel(i-td_L-td_syn),...
                W_IIret,W_IE_Ret_s,W_IE_Ret_m,W_IE_Ret_d,W_II_Ret_INs,W_IE_Ret_Rel,...
                0*I_ps(5,1,i-td_wL-td_syn),0*I_ps(5,2,i-td_wL-td_syn),kisiIret(:,i)+pnIret(:,i),zetaIret(:,i),Idc_Ret,fidN*I_dbs(2,i),n_conn_N,dt);

            rIret=rI; xIret=xI; IsIret=IsI;

    %Thalamo-cortical Relay (TCR) cells 
[vRel(:,i+1),uRel(:,i+1),rE,xE,IsE,EPSC_rel(i+1),rd,xd,Isd,EPSC_relD(i+1)]=...
    thm_rel(aErel,bErel,cErel,dErel,nErel,vRel(:,i),uRel(:,i),rErel,xErel,IsErel,rErelD,xErelD,IsErelD,...
                EPSC_rel(i-td_wL-td_syn),EPSCs(i-td_CT-td_syn),EPSCm(i-td_CT-td_syn),EPSCdF(i-td_CT-td_syn),IPSC_INs(i-td_CT-td_syn),IPSC_ret(i-td_L-td_syn),...
                W_EErel,W_EERels,W_EERelm,W_EEReld,W_EIRelINs,W_EIRelRet,...
                0*I_ps(6,1,i-td_wL-td_syn),0*I_ps(6,2,i-td_wL-td_syn),kisiErel(:,i)+pnErel(:,i),zetaErel(:,i),Idc_Rel,fidR*I_dbs(2,i),n_conn_R,dt);
            
            rErel=rE; xErel=xE; IsErel=IsE;     rErelD=rd; xErelD=xd; IsErelD=Isd;

    %Cortical Layer S
[vEs(:,i+1),uEs(:,i+1),rE,xE,IsE,EPSCs(i+1)]=...
    S_layer(aEs,bEs,cEs,dEs,nEs,vEs(:,i),uEs(:,i),rEs,xEs,IsEs,...
                EPSCs(i-td_wL-td_syn),EPSCd(i-td_L-td_syn),EPSCm(i-td_L-td_syn),EPSC_rel(i-td_TC-td_syn),IPSC_INs(i-td_wL-td_syn),IPSC_ret(i-td_TC-td_syn),...
                W_EEs,W_EEsd,W_EEsm,W_EEsRel,W_EIsINs,W_EIsRet,...
                I_ps(1,1,i-td_wL-td_syn),I_ps(1,2,i-td_wL-td_syn),kisiSE(:,i)+pnSE(:,i),zetaSE(:,i),IdcS_E,fidS*I_dbs(2,i),n_conn_S,dt);

            rEs=rE; xEs=xE; IsEs=IsE;

    %Cortical Layer M
[vEm(:,i+1),uEm(:,i+1),rE,xE,IsE,EPSCm(i+1)]=...
    M_layer(aEm,bEm,cEm,dEm,nEm,vEm(:,i),uEm(:,i),rEm,xEm,IsEm,...
                EPSCm(i-td_wL-td_syn),EPSCs(i-td_L-td_syn),EPSCd(i-td_L-td_syn),EPSC_rel(i-td_TC-td_syn),IPSC_INs(i-td_wL-td_syn),IPSC_ret(i-td_TC-td_syn),...
                W_EEm,W_EEms,W_EEmd,W_EEmRel,W_EImINs,W_EImRet,...
                I_ps(2,1,i-td_wL-td_syn),I_ps(2,2,i-td_wL-td_syn),kisiME(:,i)+pnME(:,i),zetaME(:,i),IdcM_E,fidM*I_dbs(2,i),n_conn_M,dt);
            
            rEm=rE; xEm=xE; IsEm=IsE;

    %Cortical Layer D 
[vEd(:,i+1),uEd(:,i+1),rE,xE,IsE,EPSCd(i+1),rf,xf,Isf,EPSCdF(i+1),EPSCdD(i+1)]=...
    D_layer(aEd,bEd,cEd,dEd,nEd,n_hyp,vEd(:,i),uEd(:,i),rEd,xEd,IsEd,rEdF,xEdF,IsEdF,...
                EPSCd(i-td_wL-td_syn),EPSCs(i-td_L-td_syn),EPSCm(i-td_L-td_syn),EPSC_relD(i-td_TC-td_syn),IPSC_INs(i-td_wL-td_syn),IPSC_ret(i-td_TC-td_syn),...
                W_EEd,W_EEds,W_EEdm,W_EEdRel,W_EIdINs,W_EIdRet,...
                I_ps(3,1,i-td_wL-td_syn),I_ps(3,2,i-td_wL-td_syn),kisiDE(:,i)+pnDE(:,i),zetaDE(:,i),IdcD_E,fidD*I_dbs(:,i),dt);

            rEd=rE; xEd=xE; IsEd=IsE; rEdF=rf; xEdF=xf; IsEdF=Isf;

    %Cortical Inhibitory (CI) neurons
[vIs(:,i+1),uIs(:,i+1),rI,xI,IsI,IPSC_INs(i+1)]=...
    ctx_INs(aIs,bIs,cIs,dIs,nINs,vIs(:,i),uIs(:,i),rINs,xINs,IsINs,... 
                IPSC_INs(i-td_wL-td_syn),EPSCs(i-td_wL-td_syn),EPSCm(i-td_wL-td_syn),EPSCd(i-td_wL-td_syn),EPSC_rel(i-td_TC-td_syn),IPSC_ret(i-td_TC-td_syn),...
                W_IE_INs_d,W_IE_INs_s,W_IE_INs_m,W_II_INs_Ret,W_IIins,W_IE_INs_Rel,...
                I_ps(1,1,i-td_wL-td_syn),I_ps(1,2,i-td_wL-td_syn),kisiSI(:,i)+pnSI(:,i),zetaSI(:,i),Idc_INs,fidCI*I_dbs(2,i),n_conn_CI,dt);
            
            rINs=rI; xINs=xI; IsINs=IsI;
end
ptoc;

%trash the first second of simulated APs and EPSCs (IPSCs):
vEs=vEs(:,chop_till+1:nSim); vIs=vIs(:,chop_till+1:nSim); 
vEm=vEm(:,chop_till+1:nSim); vRet=vRet(:,chop_till+1:nSim);
vEd=vEd(:,chop_till+1:nSim); vRel=vRel(:,chop_till+1:nSim);
EPSCs=EPSCs(1,chop_till+1:nSim);   IPSC_INs=IPSC_INs(1,chop_till+1:nSim);
EPSCm=EPSCm(1,chop_till+1:nSim);   IPSC_ret=IPSC_ret(1,chop_till+1:nSim);
EPSCd=EPSCd(1,chop_till+1:nSim);   EPSC_rel=EPSC_rel(1,chop_till+1:nSim);
EPSCdF=EPSCdF(1,chop_till+1:nSim); EPSC_relD=EPSC_relD(1,chop_till+1:nSim);
EPSCdD=EPSCdD(1,chop_till+1:nSim);

%deifne the PD and PD+DBS ranges
pd_ini=1; pd_fin=(nSim-chop_till)/dev;
% pd_ini=chop_till+1; pd_fin=((nSim-chop_till)/dev)+chop_till;
pd_range = pd_ini:pd_fin;
%for PD
vEsp=vEs(:,pd_range); vIsp=vIs(:,pd_range); 
vEmp=vEm(:,pd_range); vRetp=vRet(:,pd_range);
vEdp=vEd(:,pd_range); vRelp=vRel(:,pd_range);

if fidD~=0
% pd_DBS_range=pd_fin+1:(dev-1)*pd_fin;
pd_DBS_range=pd_fin+1:pd_fin+length(pd_range);
%for PD+DBS
vEsX=vEs(:,pd_DBS_range); vIsX=vIs(:,pd_DBS_range); 
vEmX=vEm(:,pd_DBS_range); vRetX=vRet(:,pd_DBS_range);
vEdX=vEd(:,pd_DBS_range); vRelX=vRel(:,pd_DBS_range);
else
pd_DBS_range=0;
vEsX=0; vIsX=0; 
vEmX=0; vRetX=0;
vEdX=0; vRelX=0;
end

clearvars -except dd Nt ij fidD fidS fidM fidCI pd_range pd_DBS_range... 
                  vEs vIs vEm vRet vEd vRel...
                  vEsp vIsp vEmp vRetp vEdp vRelp...
                  vEsX vIsX vEmX vRetX vEdX vRelX...
                  nEs nEm nEd nINs nErel nIret n_tot...
                  dt nSim tVec chop_till...
                  EPSCs EPSCm EPSCd EPSC_rel IPSC_ret IPSC_INs...
                  cIret cErel cIs cEd cEm cEs dev LFP X ECoG...
                  T dt Fs nSim td_L td_wL td_TC td_CT td_syn tVec dev n_hyp
if ij==Nt              
ptic('\n*** save the action potentials...\n ');
save(fullfile(dd,'TCM_PDwDBS_APs.mat'));
ptoc;
end
%%%%%%%%~~~~~~~~~~~~~ CONSTRUCT SPIKES from APs ~~~~~~~~~~~~~~~%%%%%%%%%%%%
if ij==Nt
ptic('\n*** evaluate spikes from the APs...\n ');
spikes=spikes_cells_wo_PS(nEs,nEm,nEd,nINs,nErel,nIret,n_tot,vEs,vEm,vEd,vIs,vRet,vRel,cIret,cErel,cIs,cEd,cEm,cEs);
if fidD~=0
spikesP=spikes_cells_wo_PS(nEs,nEm,nEd,nINs,nErel,nIret,n_tot,vEsp,vEmp,vEdp,vIsp,vRetp,vRelp,cIret,cErel,cIs,cEd,cEm,cEs);
spikesX=spikes_cells_wo_PS(nEs,nEm,nEd,nINs,nErel,nIret,n_tot,vEsX,vEmX,vEdX,vIsX,vRetX,vRelX,cIret,cErel,cIs,cEd,cEm,cEs);
end
ptoc;
end
%%%%%%%%~~~~~~~~~~~~~ CONSTRUCT LFP and ECoG signals ~~~~~~~~~~~~~~~%%%%%%%
ptic('\n*** compute LFP and ECoG signals...\n ');
% Define LFP and ECoG signals
ECoG_ij=.6*EPSCs+.2*EPSCm+.2*EPSCd-1*IPSC_INs;  
% ECoG_ij=.6*(EPSCs)+.2*(EPSCm)+.2*(EPSCd);  
LFP_ij=EPSCd-(1.0*IPSC_INs);                  %LFP definition
% LFPs=EPSCs-IPSC_INs;                %LFP only S layer
LFPd=EPSCd;                          %LFP only D layer
LFPs=EPSCs;                          %LFP only S layer
LFPm=EPSCm;                          %LFP only M layer
LFPci=-1*IPSC_INs;                          %LFP only M layer
LFPtr=-1*IPSC_ret;                          %LFP only M layer
LFPtc=EPSC_rel;                          %LFP only M layer

ECoG(ij,:)=ECoG_ij;
LFP(ij,:)=LFP_ij; 

X(1,:,ij)=LFPs; X(2,:,ij)=LFPm; X(3,:,ij)=LFPd;
X(4,:,ij)=LFPci;X(5,:,ij)=LFPtr;X(6,:,ij)=LFPtc;
ptoc;

clear vEs vIs vEm vRet vEd vRel...
      vEsp vIsp vEmp vRetp vEdp vRelp...
      vEsX vIsX vEmX vRetX vEdX vRelX...
      EPSCs EPSCm EPSCd EPSC_rel IPSC_ret IPSC_INs cIret cErel cIs cEd cEm cEs

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~CONSTRUCT RASTERS~~~~~~~~~~~~~~~~~~~~~~~~~~~
if ij==Nt
ptic('\n*** construct and save the raster plot...\n ');

%Make all rasters together in one plot (Fig. 2 and 3 in IEEE conf. paper)
Fig_spikes_all_netwrok(nEs,nEm,nEd,nINs,nErel,nIret,n_tot,spikes,dt,nSim);
filename='RasterPlot';
figdestt=fullfile(dd,filename);
saveas(gca,figdestt,'jpeg');
% savefig(figdestt); 
% saveas(gca,figdestt,'svg')
% close
ptoc;

%Make raster of the first 10 RS PNs in D
ptic('\n*** construct and save RS D raster plot...\n ');
Fig_spikes_D_RS(nIret,nErel,nINs,spikes);  
filename='D Raster RS';
figdestt=fullfile(dd,filename);
saveas(gca,figdestt,'jpeg');
savefig(figdestt); 
% saveas(gca,figdestt,'svg')
% close
ptoc;

%Make raster of the last 10 IB PNs in D
ptic('\n*** construct and save IB D raster plot...\n ');
Fig_spikes_D_IB(nIret,nErel,nINs,nEd,spikes);  
filename='D Raster IB';
figdestt=fullfile(dd,filename);
saveas(gca,figdestt,'jpeg');
savefig(figdestt); 
% saveas(gca,figdestt,'svg')
% close
ptoc;

end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if ij==Nt
ptic('\n*** save spikes and lfp signals...\n ');
save(fullfile(dd,'TCM_PDwDBS_spikes_lfp.mat'));
ptoc;
end
clearvars -except dd Nt X LFP ECoG fidD dev n_hyp...
                  T dt Fs nSim chop_till td_L td_wL td_TC td_CT td_syn tVec

% end
end

clearvars -except dd dev n_hyp fidD

%%%%%%%~~~~~~~~~~~~~` RUN the ANALYSIS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%%%%%%%
% ptic('\n*** Performing analysis: \n ');
% if fidD==0
%     Analysis_PD
% else
%     Analysis_PD_DBS
% end
% ptoc;

% clearvars
beep