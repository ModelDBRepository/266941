function [Exc_AP,Exc_Aux,r,x,Is,EPSC,rf,xf,Isf,EPSCf,EPSCd]= D_layer(a,b,c,d,n,n_hyp,v,u,r,x,Is,rf,xf,Isf,...
                                              EPSC,EPSCs,EPSCm,EPSC_rel,IPSC_in,IPSC_ret,...
                                              W_EEdd,W_EEds,W_EEdm,W_Erel,W_EI,W_Eret,... 
                                              I_psE,I_psI,kisi,zeta,Idc,Idbs,dt)
vp=30; sp=zeros(3,1);
for k=1:n
 if n_hyp==0
    Idbss=0;
 else
    if (k>=1) && (k<=n_hyp)
        Idbss=Idbs(1);
    else
        Idbss=1*Idbs(2);
    end
 end    
        v(k)=v(k)+dt*(.04*v(k)^2 + 5*v(k) -u(k) +140 +Idc(k)+... 
                   (W_EEdd(k)*EPSC/n)+(W_EEds(k)*EPSCs/n)+(W_EEdm(k)*EPSCm/n)+ ...
                   (W_EI(k)*IPSC_in/n) + (W_Eret(k)*IPSC_ret/n)+(W_Erel(k)*EPSC_rel/n) + ...    
                   I_psE - I_psI + Idbss  + kisi(k) ); 
        u(k)=u(k)+dt*a(k)*(b(k)*v(k)-u(k));
       if v(k)>= vp +zeta(k)            
%             v(k)=vp + zeta(k);   
            v(k)=c(k);               
            u(k)=u(k)+d(k);   
            sp(:)=1;
       end 
            spikeE=sp;  rr=r; xx=x; Iss=Is;
            [rs,xs,Isyn,Ipost]=TMsynE_inst(r,x,Is,spikeE);
            r=rs; x=xs; Is=Isyn;
            Ise(k)=Ipost;
            
            [rsF,xsF,IsynF,IpostF]=TMsynE_inst_F(rr,xx,Iss,spikeE);
            rf=rsF; xf=xsF; Isf=IsynF;
            IseF(k)=IpostF;
            
            [rsD,xsD,IsynD,IpostD]=TMsynE_inst_D(rr,xx,Iss,spikeE);
            rd=rsD; xd=xsD; Isd=IsynD;
            IseD(k)=IpostD;
sp(:)=0;
end
EPSC=sum(Ise); EPSCf=sum(IseF); EPSCd=sum(IseD); 
Exc_AP=v; Exc_Aux=u;