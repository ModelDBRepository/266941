 function [Exc_AP,Exc_Aux,r,x,Is,IPSC]= ctx_INs(a,b,c,d,n,v,u,r,x,Is,...
                                              IPSC,EPSCs,EPSCm,EPSCd,EPSC_rel,IPSC_ret,...
                                              W_IEd,W_IEs,W_IEm,W_Iret,W_II,W_Irel,... 
                                              I_psE,I_psI,kisi,zeta,Idc,Idbs,n_affected,dt)
vp=30; sp=zeros(3,1);
for k=1:n
 
    if (k>=1) && (k<=n_affected)
        Idbss=Idbs;
    else
        Idbss=0;
    end
        v(k)=v(k)+dt*(.04*v(k)^2 + 5*v(k) -u(k) +140 +Idc(k)+... 
                   (W_IEd(k)*EPSCd/n)+(W_IEs(k)*EPSCs/n)+(W_IEm(k)*EPSCm/n)+ ...
                   (W_II(k)*IPSC/n) + (W_Irel(k)*EPSC_rel/n)+(W_Iret(k)*IPSC_ret/n) + ...    
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
            Ise(k)=Ipost;
sp(:)=0;
end
IPSC=sum(Ise);
Exc_AP=v; Exc_Aux=u;