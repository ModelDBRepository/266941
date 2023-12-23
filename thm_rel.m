function [Exc_AP,Exc_Aux,r,x,Is,EPSC,rd,xd,Isd,EPSC_d]= thm_rel(a,b,c,d,n,v,u,r,x,Is,rd,xd,Isd,...
                                                   EPSC,EPSCs,EPSCm,EPSCd,IPSC_in,IPSC_ret,...
                                                   W_EE,W_EErs,W_EErm,W_EErd,W_EI_IN,W_EI_ret,...
                                                   I_psE,I_psI,kisi,zeta,Idc,Idbs,n_affected,dt)
vp=30; sp=zeros(3,1);
for k=1:n
 
    if (k>=1) && (k<=n_affected)
        Idbss=Idbs;
    else
        Idbss=0;
    end
    v(k)=v(k)+dt*(.04*v(k)^2 + 5*v(k) -u(k) +140 +Idc(k)+...
         (W_EE(k)*EPSC/n)+ ...                                                 %self feedback
         (W_EErs(k)*EPSCs/n) + (W_EErm(k)*EPSCm/n) + (W_EErd(k)*EPSCd/n) +... %Excitatory inputs from S, M and D layers
         (W_EI_ret(k)*IPSC_ret/n) + (W_EI_IN(k)*IPSC_in/n)+ ...                %Inhibitory inputs from INs and reticulars
          I_psE - I_psI + Idbss + kisi(k) ); 
    u(k)=u(k)+dt*a(k)*(b(k)*v(k)-u(k));
      if v(k)>= vp +zeta(k)            
           v(k)=vp + zeta(k);   
           v(k)=c(k);             
           u(k)=u(k)+d(k);   
           sp(:)=1;
      end 
  spikeE=sp; rr=r; xx=x; Iss=Is;
      [rs,xs,Isyn,Ipost]=TMsynE_inst(r,x,Is,spikeE);
      r=rs; x=xs; Is=Isyn;
      Ise(k)=Ipost;
  
      [rsD,xsD,IsynD,IpostD]=TMsynE_inst_D(rr,xx,Iss,spikeE);
      rd=rsD; xd=xsD; Isd=IsynD;
      IseD(k)=IpostD;
sp(:)=0;
end

EPSC=sum(Ise); EPSC_d=sum(IseD);
Exc_AP=v; Exc_Aux=u;