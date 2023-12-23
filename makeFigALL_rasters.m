function fig=makeFigALL_rasters(nEs,nEm,nEd,nIs,nIret,nErel,spike_times_Es,spike_times_Em,spike_times_Ed,spike_times_Is,spike_times_Iret,spike_times_Erel,vEs,vEm,vEd,vIs,vIret,vErel,xaxis,dt,nSim)

fig=figure('units','normalized','outerposition',[0 0 1 1]) ;
% fig.OuterPosition=[230 250 570 510];

ax1=subplot(6,2,1); 
    title('Cortical Excitatory Neurons in Layer S'); hold on
for ik=1:nEs
for ij=1:length(spike_times_Es{ik,:})
    line([spike_times_Es{ik}(ij) spike_times_Es{ik}(ij)],[ik-1 ik],'LineWidth',0.5)%, 'Color','k','LineWidth',1.5)
end
% axis([0 nSim 0 nEs]); hold on;
end
ylim([0 nEs])
ylabel('E neurons'); hold on
set(gca,'FontSize',12,'FontWeight','bold')

ax2=subplot(6,2,2); title('Action potential examples'); hold on 
plot(dt*xaxis,squeeze(vEs(1,xaxis)),'LineWidth',0.5); hold on; ylabel('v_E1s(mV)')
ylim([-100 100])
set(gca,'FontSize',12,'FontWeight','bold')

ax3=subplot(6,2,3); 
    title('Cortical Excitatory Neurons in Layer M'); hold on
for ik=1:nEm
for ij=1:length(spike_times_Em{ik,:})
    line([spike_times_Em{ik}(ij) spike_times_Em{ik}(ij)],[ik-1 ik],'LineWidth',0.5)%, 'Color','k','LineWidth',1.5)
end
% axis([0 nSim 0 nEm]); hold on;
end
ylim([0 nEs])
ylabel('E neurons'); hold on
set(gca,'FontSize',12,'FontWeight','bold')

ax4=subplot(6,2,4);
plot(dt*xaxis,squeeze(vEm(1,xaxis)),'LineWidth',0.5); hold on; ylabel('v_E1m(mV)')
ylim([-100 100])
set(gca,'FontSize',12,'FontWeight','bold')

ax5=subplot(6,2,5); 
    title('Cortical Excitatory Neurons in Layer D'); hold on
for ik=1:nEd
for ij=1:length(spike_times_Ed{ik,:})
    line([spike_times_Ed{ik}(ij) spike_times_Ed{ik}(ij)],[ik-1 ik],'LineWidth',0.5)%, 'Color','k','LineWidth',1.5)
end
% axis([0 nSim 0 nEd]); hold on;
end
ylim([0 nEs])
ylabel('E neurons'); hold on
set(gca,'FontSize',12,'FontWeight','bold')

ax6=subplot(6,2,6);
plot(dt*xaxis,squeeze(vEd(1,xaxis)),'LineWidth',0.5); hold on; ylabel('v_E1d(mV)')
ylim([-100 100])
set(gca,'FontSize',12,'FontWeight','bold')


ax7=subplot(6,2,7);
title('Cortical Inhibitory Interneurons'); hold on
for ik=1:nIs
for ij=1:length(spike_times_Is{ik,:})
    line([spike_times_Is{ik}(ij) spike_times_Is{ik}(ij)],[ik-1 ik],'LineWidth',0.5)
end
% axis([0 nSim 0 nIs]); hold on;
end
ylim([0 nEs])
ylabel('I neurons'); 
set(gca,'FontSize',12,'FontWeight','bold')

ax8=subplot(6,2,8);
plot(dt*xaxis,squeeze(vIs(1,xaxis)),'LineWidth',0.5); ylabel('v_I1in(mV)'); zoom xon
ylim([-100 100])
set(gca,'FontSize',12,'FontWeight','bold')

ax9=subplot(6,2,9);
title('Thalamic Reticular Neurons'); hold on
for ik=1:nIret
for ij=1:length(spike_times_Iret{ik,:})
    line([spike_times_Iret{ik}(ij) spike_times_Iret{ik}(ij)],[ik-1 ik],'LineWidth',0.5)
end
% axis([0 nSim 0 nIret]); hold on;
end
ylim([0 nEs])
ylabel('I neurons');  
set(gca,'FontSize',12,'FontWeight','bold')

ax10=subplot(6,2,10);
plot(dt*xaxis,squeeze(vIret(1,xaxis)),'LineWidth',0.5); ylabel('v_I1ret(mV)'); zoom xon
ylim([-100 100])
set(gca,'FontSize',12,'FontWeight','bold')


ax11=subplot(6,2,11);
title('Thalamocortical Relay Neurons'); hold on
for ik=1:nErel
for ij=1:length(spike_times_Erel{ik,:})
    line([spike_times_Erel{ik}(ij) spike_times_Erel{ik}(ij)],[ik-1 ik],'LineWidth',0.5)
end
% axis([0 nSim 0 nErel]); hold on;
end
ylim([0 nEs])
xlabel('Time (s)');
ylabel('E neurons');  
set(gca,'FontSize',12,'FontWeight','bold')
zoom xon;

ax12=subplot(6,2,12);
plot(dt*xaxis,squeeze(vErel(1,xaxis)),'LineWidth',0.5); ylabel('v_E1rel(mV)'); zoom xon
linkaxes([ax1,ax3,ax5,ax7,ax9,ax11],'x'); 
xlim([0 dt*nSim]); ylim([-100 100])
linkaxes([ax2,ax4,ax6,ax8,ax10,ax12],'x'); xlabel('Time (ms)'); 
xlim([0 dt*nSim]); ylim([-100 100])
set(gca,'FontSize',12,'FontWeight','bold')