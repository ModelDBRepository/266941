function fig_DSpikes=Fig_spikes_D_RS(n_ret,n_rel,n_CI,spike_times)  

fig_DSpikes=figure; set(gcf,'Visible','off');
title('10 RS PNs in D'); hold on

for ik=n_ret+n_rel+n_CI+1:n_ret+n_rel+n_CI+10
for ij=1:length(spike_times{ik,:})
    line([spike_times{ik}(ij) spike_times{ik}(ij)],[ik-1 ik],'LineWidth',1,'Color','k')
end
% axis([0 nSim 0 n_t]); hold on;
end

% xlim([0 dt*nSim]); zoom xon
% ylim([0 n_t])

xlabel('Time (ms)');
ylabel('Neuron #'); hold on
set(gca,'FontSize',12,'FontWeight','bold')