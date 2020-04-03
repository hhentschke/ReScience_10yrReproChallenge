% plots histogram of peak amplitude and inter-peak-intervals (ipi) of
% theta. needs exported data from rmouse_oscillRegularity.m

% **** plots are intended to be combined with those of rawnstreams03 - make
% sure the important settings (axis limits, reduction factors, etc.) are
% identical ****

% *** Sep 2010: code will not work, names of variables exported by
% rmouse_oscillRegularity have changed (but should run again with minor
% adjustments)

ornt='portrait';
labelscale('fontSz',8,'scaleFac',.5,'lineW',.75,'markSz',4);

figdir='d:\projects\rmouse\paper_atropine\rawFig\';
% figdir='';
printas=[];
% printas='-dpsc2';
figName=mfilename;
figure(1), clf, orient(ornt)

% make sure this is the current y lim of rawnstreams03
yl=[-1.5 1.5];

rmouse_ini;

fn={'04708001_proc_exc1_oscillReg.mat','04708002_proc_exc1_oscillReg.mat'};

for fi=1:length(fn)
  load([figdir fn{fi}]);
  chaI=12;

  % immobile first
  sph(1)=subplot(2,3,1); hold on
  % ph=stairs(binA,tmpr(3).thNegPeakAH(:,chaI),'k');
  ph=contourbarh(binA,tmpr(3).thNegPeakAH(:,chaI),'color','k');
  if fi==1
    set(ph,'color',[.6 .6 .6]);
  end
  axis tight
  
  % exploring
  sph(2)=subplot(2,3,4); hold on
  ph=contourbarh(binA,tmpr(4).thNegPeakAH(:,chaI),'color','k');
  if fi==1
    set(ph,'color',[.6 .6 .6]);
  end

end
subpax(gcf);
set(sph,'ylim',yl);



if ~isempty(printas), 
  print(printas,[figdir figName]); 
end
