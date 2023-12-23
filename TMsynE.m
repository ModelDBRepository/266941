function [Ipost,tVec]=TMsynE(tevent,nSim,tdelay,dt)

tVec = tdelay+1:nSim-1;

tauf=[670,17,326];
taud=[138,671,329];
U=[.09,.5,.29];
A=[.20,.63,.17].*ones(1,3);
% A=[.63,.30,.07].*ones(1,3);
% A=[.5,.5,0].*ones(1,3);

%I Kernel time constant
% tausE=1.75; %*
tausE=3;

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
       Is(p,i+1) = Is(p,i) + dt*((-1/tausE)*Is(p,i) + A(p)*r(p,i+1)*x(p,i)*(spd(i-tdelay)));
   end
end
    
Ipost=sum(Is,1);