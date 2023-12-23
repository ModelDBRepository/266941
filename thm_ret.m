function [Inh_AP,Inh_Aux,r,x,Is,IPSC]= thm_ret(a,b,c,d,n,v,u,r,x,Is,...
                                                   IPSC,EPSCs,EPSCm,EPSCd,IPSC_in,EPSC_rel,...
                                                   W_II,W_IErs,W_IErm,W_IErd,W_II_IN,W_IE_rel,...
                                                   I_psE,I_psI,kisi,zeta,Idc,Idbs,n_affected,dt)
vp=30; sp=zeros(3,1);       
for k=1:n
    if (k>=1) && (k<=n_affected)
        Idbss=Idbs;
    else
        Idbss=0;
    end
    v(k)=v(k)+dt*(.04*v(k)^2 + 5*v(k) -u(k) +140 +Idc(k)+...
         (W_II(k)*IPSC/n)+ ...                                                                          %self feedback
         (W_IErs(k)*EPSCs/n) + (W_IErm(k)*EPSCm/n) + (W_IErd(k)*EPSCd/n) +(W_IE_rel(k)*EPSC_rel/n)+ ... %Excitatory inputs from S, M and D layers plus relay cells
         (W_II_IN(k)*IPSC_in/n)+ ...                                                                    %Inhibitory inputs from INs
          I_psE - I_psI + Idbss + kisi(k)); 
    u(k)=u(k)+dt*a(k)*(b(k)*v(k)-u(k));
      if v(k)>= vp +zeta(k)            
           v(k)=vp + zeta(k);   
           v(k)=c(k);             
           u(k)=u(k)+d(k);   
           sp(:)=1;
      end 
  spikeI=sp;
     [rs,xs,Isyn,Ipost]=TMsynI_inst(r,x,Is,spikeI);
      r=rs; x=xs; Is=Isyn;
  Isi(k)=Ipost;
sp(:)=0;
end

IPSC=sum(Isi);
Inh_AP=v; Inh_Aux=u;