function [r_new,x_new,Isyn,Ipost]=TMsynE_inst_F(r,x,Is,sp_event)
dt=.1;

tauf=[670,17,326];
taud=[138,671,329];
U=[.09,.5,.29];
A=[1,0,0].*ones(1,3);

%I Kernel time constant
% taus=1.75; %*
taus=3;

for p=1:3 %F,D, P synapses
       r(p) = r(p) + dt*(-(r(p)/tauf(p))+U(p)*(1-r(p))*sp_event(p)); 
       x(p) = x(p) + dt*((1/taud(p))*(1-x(p)) - (r(p)+U(p)*(1-r(p)))*x(p)*sp_event(p));
       Is(p) = Is(p) + dt*((-1/taus)*Is(p) + A(p)*(r(p)+U(p)*(1-r(p)))*x(p)*sp_event(p));
end

r_new=r;
x_new=x;
Isyn=Is;
Ipost=sum(Is);
