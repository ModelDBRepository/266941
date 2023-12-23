function figg=lfp_fig(LFP,LFP_fil,f0,p0,f1,p1,f2,p2,f3,p3,dt,tVec,fid)

figg=figure('units','normalized','outerposition',[0 0 1 1]) ; 
% figg.OuterPosition=[0 250 2030 300];

ax1=subplot(3,2,[1,2]); title('Raw LFP'); hold on
plot(dt*tVec(1:length(LFP)),LFP,'LineWidth',1); zoom xon; 
xlim([0 dt*length(LFP)])
set(gca,'FontSize',12,'FontWeight','bold')

ax2=subplot(3,2,[3,4]); title('LFP, \beta filterred'); hold on
plot(dt*tVec(1:length(LFP)),LFP_fil,'LineWidth',1)
xlabel('Time (ms)');
set(gca,'FontSize',12,'FontWeight','bold')
xlim([0 dt*length(LFP)])
linkaxes([ax1,ax2],'x')

ax3=subplot(3,2,5); title('Not filterred'); hold on
plot(f0,p0,'Linewidth',1); hold on 
plot(f1,p1,'Linewidth',1); zoom xon
if fid==0
xlim([0 100]);
else
xlim([0 300]);
end
xlabel('Frequency (Hz)')
% ylim([0 .001])
ylabel('Power')
set(gca,'FontSize',12,'FontWeight','bold')

ax4=subplot(3,2,6); title('\beta filterred'); hold on
plot(f2,p2,'Linewidth',1); hold on
plot(f3,p3,'Linewidth',1); zoom xon
xlim([0 100]); 
xlabel('Frequency (Hz)')
% ylim([0 .25])
%ylim([-3 1])
set(gca,'FontSize',12,'FontWeight','bold')