function fig_AllSpikes=Fig_spikes_all_netwrok(n_s,n_m,n_d,n_IN,n_rel,n_ret,n_t,spike_times,dt,nSim)  

fig_AllSpikes=figure; set(gcf,'Visible','off');
title('Thalamo-Cortical netwrok'); hold on

% obj1=patch([0 dt*nSim dt*nSim 0],[0 0 n_ret n_ret],'b'); hold on
% obj2=patch([0 dt*nSim dt*nSim 0],[n_ret n_ret n_ret+n_rel n_ret+n_rel],'r'); hold on
% obj3=patch([0 dt*nSim dt*nSim 0],[n_ret+n_rel n_ret+n_rel n_ret+n_rel+n_IN n_ret+n_rel+n_IN],'y'); hold on
% obj4=patch([0 dt*nSim dt*nSim 0],[n_ret+n_rel+n_IN n_ret+n_rel+n_IN n_ret+n_rel+n_IN+n_d n_ret+n_rel+n_IN+n_d],'g'); hold on
% obj5=patch([0 dt*nSim dt*nSim 0],[n_ret+n_rel+n_IN+n_d n_ret+n_rel+n_IN+n_d n_ret+n_rel+n_IN+n_d+n_m n_ret+n_rel+n_IN+n_d+n_m],'m'); hold on
% obj6=patch([0 dt*nSim dt*nSim 0],[n_ret+n_rel+n_IN+n_d+n_m n_ret+n_rel+n_IN+n_d+n_m n_ret+n_rel+n_IN+n_d+n_m+n_s n_ret+n_rel+n_IN+n_d+n_m+n_s],'c'); hold on
% transp=.03;
% alpha(obj1,transp); alpha(obj2,transp); alpha(obj3,transp); alpha(obj4,transp); alpha(obj5,transp); alpha(obj6,transp);

for ik=1:n_t
for ij=1:length(spike_times{ik,:})
    line([spike_times{ik}(ij) spike_times{ik}(ij)],[ik-1 ik],'LineWidth',1,'Color','k')
end
% axis([0 nSim 0 n_t]); hold on;
end

% xlim([0 dt*nSim]); zoom xon
ylim([0 n_t])

% line([0 nSim],[n_ret n_ret],'LineWidth',.5,'Color','k','LineStyle','-'); hold on
% line([0 nSim],[n_ret+n_rel n_ret+n_rel],'LineWidth',.5,'Color','k','LineStyle','-'); hold on
% line([0 nSim],[n_ret+n_rel+n_IN n_ret+n_rel+n_IN],'LineWidth',.5,'Color','k','LineStyle','-'); hold on
% line([0 nSim],[n_ret+n_rel+n_IN+n_s n_ret+n_rel+n_IN+n_s],'LineWidth',.5,'Color','k','LineStyle','-'); hold on
% line([0 nSim],[n_ret+n_rel+n_IN+n_s+n_m n_ret+n_rel+n_IN+n_s+n_m],'LineWidth',.5,'Color','k','LineStyle','-'); hold on


% dim1 = [.7 .1]; str1= 'Thalamic Reticular Nucleus (TRN)'; annotation('textbox',dim1,'String',str1,'BackgoundColor','b','HorizontalAlignment','center','VerticalAlignment','middle'); drawnow
% dim2 = [.7 .2]; str2= 'Thalamocortical Relacy Nucleus (TCR)'; annotation('textbox',dim2,'String',str2,'BackgoundColor','r','HorizontalAlignment','center','VerticalAlignment','middle'); drawnow
% dim3 = [.7 .4]; str3= 'Cortical Inhibitory Neurons (IN)'; annotation('textbox',dim3,'String',str3,'BackgoundColor','b','HorizontalAlignment','center','VerticalAlignment','middle'); drawmow
% dim4 = [.7 .6]; str4= 'Deep Cortical Layer (D)'; annotation('textbox',dim4,'String',str4,'BackgoundColor','g','HorizontalAlignment','center','VerticalAlignment','middle'); drawnow
% dim5 = [.7 .8]; str5= 'Middle Cortical Layer (M)'; annotation('textbox',dim5,'String',str5,'BackgoundColor','y','HorizontalAlignment','center','VerticalAlignment','middle'); drawnow
% dim6 = [.7 .9]; str6= 'Superficial Cortical Layer (S)'; annotation('textbox',dim6,'String',str6,'BackgoundColor','c','HorizontalAlignment','center','VerticalAlignment','middle'); drawnow


xlabel('Time (ms)');
ylabel('Neuron #'); hold on
set(gca,'FontSize',12,'FontWeight','bold')