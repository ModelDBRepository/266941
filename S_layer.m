function [Exc_AP,Exc_Aux,r,x,Is,EPSC]= S_layer(a,b,c,d,n,v,u,r,x,Is,...
                                              EPSC,EPSCd,EPSCm,EPSC_rel,IPSC_in,IPSC_ret,...
                                              W_EEss,W_EEsd,W_EEsm,W_Erel,W_EI,W_Eret,... 
                                              I_psE,I_psI,kisi,zeta,Idc,Idbs,n_affected,dt)
vp=30; sp=zeros(3,1);
for k=1:n
 
    if (k>=1) && (k<=n_affected)
        Idbss=Idbs;
    else
        Idbss=0;
    end
        v(k)=v(k)+dt*(.04*v(k)^2 + 5*v(k) -u(k) +140 +Idc(k)+... 
                   (W_EEss(k)*EPSC/n)+(W_EEsd(k)*EPSCd/n)+(W_EEsm(k)*EPSCm/n)+ ...
                   (W_EI(k)*IPSC_in/n) + (W_Eret(k)*IPSC_ret/n)+(W_Erel(k)*EPSC_rel/n) + ...    
                   I_psE - I_psI + Idbss + kisi(k) ); 
        u(k)=u(k)+dt*a(k)*(b(k)*v(k)-u(k));
       if v(k)>= vp +zeta(k)         
            v(k)=vp + zeta(k);   
            v(k)=c(k);               
            u(k)=u(k)+d(k);   
            sp(:)=1;
       end 
            spikeE=sp;
            [rs,xs,Isyn,Ipost]=TMsynE_inst(r,x,Is,spikeE);
            r=rs; x=xs; Is=Isyn;
            Ise(k)=Ipost;
sp(:)=0;
end
EPSC=sum(Ise);
Exc_AP=v; Exc_Aux=u;