% Devision Factor:
% fac=2.5;
fac=fac_N;

ini=0.0; fin=1; int=fin-ini; 

r_s=ini+int*rand(nEs,1);     %These are to restrict the normalized distribution variance or deviation from the mean
r_m=ini+int*rand(nEm,1);
r_d=ini+int*rand(nEd,1);
r_ins=ini+int*rand(nINs,1);
r_ret=ini+int*rand(nIret,1);
r_rel=ini+int*rand(nErel,1); 

%COUPLING STRENGTHs within each structure (The same in Normal and PD)
aee_s=-1e1/fac;       W_EEs= aee_s*r_s;             %Layer S (was -1e-2 for IEEE paper)
aee_m=-1e1/fac;       W_EEm= aee_m*r_m;             %Layer M (was -1e-2 for IEEE paper)
aee_d=-1e1/fac;       W_EEd= aee_d*r_d;             %Layer D (was -1e-2 for IEEE paper)
aii_INs=-5e2/fac;    W_IIins= aii_INs*r_ins;        %INs
aii_ret=-5e1/fac;    W_IIret= aii_ret*r_ret;        %Reticular cells
aee_rel=0/fac;       W_EErel= aee_rel*r_rel;        %Relay cells

%COUPLING STRENGTHs between structures (PD)
%S
aee_sm=1e1/fac;   W_EEsm=aee_sm*r_s;             %M to S couplings
aee_sd=5e2/fac;      W_EEsd=aee_sd*r_s;             %D to S couplings
aei_sINs=-5e2/fac;    W_EIsINs=aei_sINs*r_s;         %INs to S couplings
aei_sRet=0/fac;      W_EIsRet=aei_sRet*r_s;         %Ret. to S couplings
aee_sRel=0/fac;      W_EEsRel=aee_sRel*r_s;         %Rel. to S couplings
%M
aee_ms=3e2/fac;      W_EEms=aee_ms*r_m;             %S to M couplings
aee_md=0/fac;        W_EEmd=aee_md*r_m;             %D to M couplings
aei_mINs=-3e2/fac;    W_EImINs=aei_mINs*r_m;         %INs to M couplings
aei_mRet=0/fac;      W_EImRet=aei_mRet*r_m;         %Ret. to M couplings
aee_mRel=0/fac;      W_EEmRel=aee_mRel*r_m;         %Rel. to M couplings
%D
aee_ds=3e2/fac;      W_EEds=aee_ds*r_d;             %S to D couplings
aee_dm=0/fac;        W_EEdm=aee_dm*r_d;             %M to D couplings
aei_dINs=-7.5e3/fac;   W_EIdINs=aei_dINs*r_d;         %INs to D couplings
aei_dRet=0/fac;      W_EIdRet=aei_dRet*r_d;         %Ret. to D couplings
aee_dRel=1e1/fac;    W_EEdRel=aee_dRel*r_d;         %Rel. to D couplings
%INs
aie_inss=2e2/fac;    W_IE_INs_s=aie_inss*r_ins;     %S to INs couplings
aie_insm=2e2/fac;    W_IE_INs_m=aie_insm*r_ins;     %M to INs couplings
aie_insd=2e2/fac;    W_IE_INs_d=aie_insd*r_ins;     %D to INs couplings
aii_InsRet=0/fac;    W_II_INs_Ret=aii_InsRet*r_ins; %Ret. to INs couplings
aie_InsRel=1e1/fac;  W_IE_INs_Rel=aie_InsRel*r_ins; %Rel. to INs couplings
%Ret.
aie_rets=0/fac;      W_IE_Ret_s=aie_rets*r_ret;     %S to Ret couplings
aie_retm=0/fac;      W_IE_Ret_m=aie_retm*r_ret;     %M to Ret couplings
aie_retd=7e2/fac;    W_IE_Ret_d=aie_retd*r_ret;     %D to Ret couplings
aii_RetIns=0/fac;    W_II_Ret_INs=aii_RetIns*r_ret; %Ret. Ret INs couplings
aie_RetRel=1e3/fac;  W_IE_Ret_Rel=aie_RetRel*r_ret; %Rel. Ret INs couplings
%Rel.
aee_rels=0/fac;      W_EERels=aee_rels*r_rel;       %S to Rel couplings
aee_relm=0/fac;      W_EERelm=aee_relm*r_rel; 	    %M to Rel couplings
aee_reld=7e2/fac;    W_EEReld=aee_reld*r_rel;       %D to Rel couplings
aei_RelINs=0/fac;    W_EIRelINs=aei_RelINs*r_rel;   %INs to Rel couplings
aei_RelRet=-5e2/fac; W_EIRelRet=aei_RelRet*r_rel;   %Ret to Rel couplings
