% scaling et al for figures displaying lags of the cross correlations
set(axh(1),'xlim',[-20 65],'xtick',-20:20:80);
set(axh(2),'xlim',[-20 65],'xtick',-20:20:80);
subplot(axh(1));
lh=line([0 0],get(gca,'ylim'),'linestyle','--','linewidth',get(ph(1),'linewidth')*.5,'color','k');
subplot(axh(2));
lh=line([0 0],get(gca,'ylim'),'linestyle','--','linewidth',get(ph(1),'linewidth')*.5,'color','k');


set(axh(1),'xlim',[-1.2 3.5],'xtick',-1:3);
set(axh(2),'xlim',[-1.2 3.5],'xtick',-1:3);
subplot(axh(1));
lh=line([0 0],get(gca,'ylim'),'linestyle','--','linewidth',get(ph(1),'linewidth')*.5,'color','k');
subplot(axh(2));
lh=line([0 0],get(gca,'ylim'),'linestyle','--','linewidth',get(ph(1),'linewidth')*.5,'color','k');


set(axh(1),'xlim',[0 40],'xtick',0:10:40);
set(axh(2),'xlim',[0 40],'xtick',0:10:40);

