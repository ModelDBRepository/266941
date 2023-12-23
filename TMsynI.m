function [Ipost,tVec]=TMsynI(tevent,nSim,tdelay,dt)

tVec = tdelay+1:nSim-1;

tauf=[376,21,62];
taud=[45,706,144];
U=[.016,.25,.32];
A=[.08,.75,.17].*ones(1,3);
%I Kernel time constant
tausI=11;  %*
% tausI=3;

% Initial values
r=zeros(3,nSim);
x=ones(3,nSim);
Is=zeros(3,nSim);
spd=zeros(1,nSim); 

for p=1:3 %F,D, P synapses
spd(tevent)=1/dt;   %for poissonian input
   for i=tdelay+1:nSim-1  
       r(p,i+1) = r(p,i) + dt*(-(r(p,i)/tauf(p))+U(p)*(1-r(p,i))*spd(i-tdelay)); 
       x(p,i+1) = x(p,i) + dt*((1/taud(p))*(1-x(p,i)) - r(p,i+1)*x(p,i)*spd(i-tdelay));
       Is(p,i+1) = Is(p,i) + dt*((-1/tausI)*Is(p,i) + A(p)*r(p,i+1)*x(p,i)*(spd(i-tdelay)));
   end
end
    
Ipost=sum(Is,1);