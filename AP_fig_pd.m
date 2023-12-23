function figx=AP_fig_pd(v_S,v_M,v_D,v_CI,v_TCR,v_TRN)

figx=figure;
title('Action potential examples'); hold on

subplot(6,1,1)
% title('S'); hold on
for nc=1:10
plot((nc-1)*200+v_S(nc,:),'k');hold on
end
% xlim([10000 20000]); 
ylim([-200 nc*190])
axis off

subplot(6,1,2)
% title('M'); hold on
for nc=1:10
plot((nc-1)*200+v_M(nc,:),'k');hold on
end
% xlim([10000 20000]); 
ylim([-200 nc*190])
axis off

subplot(6,1,3)
% title('D'); hold on
for nc=1:10
plot((nc-1)*200+v_D(nc,:),'k');hold on
end
% xlim([10000 20000]); 
ylim([-200 nc*190])
axis off

subplot(6,1,4)
% title('CI','Position','left'); hold on
for nc=1:10
plot((nc-1)*200+v_CI(nc,:),'k');hold on
end
% xlim([10000 20000]); 
ylim([-200 nc*190])
axis off

subplot(6,1,5)
% title('TCR'); hold on
for nc=1:10
plot((nc-1)*200+v_TCR(nc,:),'k');hold on
end
% xlim([10000 20000]); 
ylim([-200 nc*190])
axis off

subplot(6,1,6)
% title('TRN'); hold on
for nc=1:4
plot((nc-1)*200+v_TRN(nc,:),'k');hold on
end
% xlim([10000 20000]); 
ylim([-200 nc*190])
axis off