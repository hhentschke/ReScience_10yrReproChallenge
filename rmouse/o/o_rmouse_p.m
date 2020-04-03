% generate contour plots of CC?
cc_contourP=1;
% generate contour plots spec analysis?
spec_contourP=1;

% which 'segmentTypes (=behavior classes) to plot?
bpix=[immV explV];

% some shorties
pcIdx=AP.pcIdx;
pcInd=AP.pcInd;
LFPpcInd1=AP.LFPpcInd1;
LFPpcInd2=AP.LFPpcInd2;
LFPccInd=AP.LFPccInd;
% template
tempo=repmat(nan,nAllLFPCh,1);
% color associated with principal channel
pcCol=[1 .6 .6];
% load r only if it does not yet exist
if ~exist('r','var'),
  if exist([AP.resFn '.mat'],'file')
    load(AP.resFn);
  else
    logstr{end+1}='results file does not exist';
    warning(logstr{end});
    r=[];
  end
end
% damn..
if strncmp(WP.mver,'6.',2)
  cbPos='vert';
else
  cbPos='EastOutside';
end

% warnings like 'negative data ignored' etc keep popping up, so kill them
warning off


% -------------------------------------------------------------------------------
%  ***  FIGURE I: general data overview & summary plots of spectral analysis ***
% -------------------------------------------------------------------------------
figure(sumFigH);
% plots showing spectral data:
% ------ power spectral densities (no cross sd) from raw data:
if isfield(r,'rawPMn')
  % -- exploring, all electrodes
  subplot('position',[.33-.33/2+xmarg .8/3*2+ymarg .33/2-2*xmarg .8/3-2*ymarg]);
  hold on
  if length(diag(r(explV).rawPMn))>1
    % plot from 1 to 180 Hz
    frix=r(1).F>=1 & r(1).F<=180;
    % extract spectra, omitting non analyzed channels (nans)
    tmpP=cat(2,r(explV).rawPMn{AP.dixie(LFPccInd)});
    plot(r(1).F(frix),tmpP(frix,:),'k');
    hold on;
    % plot principal channel in red
    ph=plot(r(1).F(frix),tmpP(frix,LFPpcInd1),'k'); % sic! LFPpcInd1
    set(ph,'color',pcCol);
    % for whatever reason loglog doesnt work, so doit by hand
    set(gca,'xscale','log','yscale','log');
    % niceyax etc does not work with log scaling
    axis tight
    yl=get(gca,'ylim');
    % set(gca,'ylim',[10^(floor(log10(yl(1)))) 10^(ceil(log10(yl(2))))])
    ylabel('power spectral density (mV^2/Hz)');
    xlabel('frequency (Hz)');
    th=title(r(explV).segmentType); set(th,'fontsize',10,'fontweight','bold');
  end
end

if isfield(r,'rawPMn')
  % -- principal electrode, different behaviors
  subplot('position',[.33-.33/2+xmarg .8/3+ymarg .33/2-2*xmarg .8/3-2*ymarg]);
  hold on
  for i=bpix
    if length(diag(r(i).rawPMn))>1
      pcol=AP.segmentType{i,3};
      frix=r(1).F>=1 & r(1).F<=180;
      ph=plot(r(1).F(frix),r(i).rawPMn{LFPpcInd2,LFPpcInd2}(frix));
      set(ph,'color',pcol);
      hold on;
      set(gca,'xscale','log','yscale','log','xtick',[10 100]);
    end
    % niceyax etc does not work with log scaling
    axis tight
    ylabel('power spectral density (mV^2/Hz)');
    xlabel('frequency (Hz)');
    th=title('princ channel'); set(th,'fontsize',10,'fontweight','bold');
  end
end

% ------ power in freq bands (loop over bands):
for sigi=1:4
  switch sigi
    case 1
      fifi='rawDePEMn';
      ppos=[.33+xmarg .8/4*3+ymarg .33/2-2*xmarg .8/4-2*ymarg];
      titulo='{\delta}, Power ';
    case 2
      fifi='rawThPEMn';
      ppos=[.33+xmarg .8/4*2+ymarg .33/2-2*xmarg .8/4-2*ymarg];
      titulo='{\theta}_{narrow}, Power ';
    case 3
      fifi='rawGaPEMn';
      ppos=[.33+xmarg .8/4*1+ymarg .33/2-2*xmarg .8/4-2*ymarg];
      titulo='{\gamma}, Power';
    case 4
      fifi='rawRiPEMn';
      ppos=[.33+xmarg .8/4*0+ymarg .33/2-2*xmarg .8/4-2*ymarg];
      titulo='Rip, Power';
  end
  ylab='P (mV^2)';
  if isfield(r,fifi)
    subplot('position',ppos);
    hold on
    for i=bpix
      % assume that if theta power exists power had also been calculated in other bands
      if length(diag(r(i).rawThNarrowPEMn))>1
        % more handy constants..
        pcol=AP.segmentType{i,3}; psym=AP.segmentType{i,4};
        if sigi==1
          tmpMn=cat(2,r(i).rawDePEMn{AP.dixie});
          tmpStd=cat(2,r(i).rawDePEStd{AP.dixie});
        elseif sigi==2
          tmpMn=cat(2,r(i).rawThNarrowPEMn{AP.dixie});
          tmpStd=cat(2,r(i).rawThNarrowPEStd{AP.dixie});
        elseif sigi==3
          tmpMn=cat(2,r(i).rawGaPEMn{AP.dixie});
          tmpStd=cat(2,r(i).rawGaPEStd{AP.dixie});
        else
          tmpMn=cat(2,r(i).rawRiPEMn{AP.dixie});
          tmpStd=cat(2,r(i).rawRiPEStd{AP.dixie});
        end
        ph=errorbar(WP.elx,tmpMn,zeros(size(tmpMn)),tmpStd,[psym '-']);
        set(ph,'color',pcol);
      end
    end
    % line (below) will not be drawn if lower y limit is <0 and then yscale set to
    % log, so set y axis limits manually
    nicexyax;
    % yl=get(gca,'ylim');
    % yl=yl.*[.8 1.2];
    % set(gca,'ylim',yl,'yscale','log','xtick',WP.elx,'xticklabel',WP.xtl1);
    set(gca,'xtick',WP.elx,'xticklabel',WP.xtl1);
    ylabel(ylab);
    % line to indicate principal electrode
    lh=line(WP.elx(LFPpcInd2)*[1 1]',get(gca,'ylim'));
    set(lh,'color',pcCol,'linestyle',':');
    % push to background
    set(gca,'children',circshift(get(gca,'children'),-1))
    % second axis for electrode #
    ax2=axes('position',get(gca,'position'),'color','none','xaxisloc','top',...
      'yaxisloc','right','ylim',get(gca,'ylim'),'ytick',get(gca,'ytick'),...
      'xlim',get(gca,'xlim'),'xtick',WP.elx,'xticklabel',WP.xtl2,'yticklabel',[]);
    th=title(titulo); set(th,'fontsize',10,'fontweight','bold');
    tp=get(th,'position');
    yl=get(gca,'ylim');
    tp(2)=yl(1)+diff(yl)*1.085;
    set(th,'position',tp);
  end
end

% ------ coherence in freq bands (loop over bands):
for sigi=1:4
  switch sigi
    case 1
      fifi='rawCohMnDe';
      ppos=[.33+.33/2+xmarg .8/4*3+ymarg .33/2-2*xmarg .8/4-2*ymarg];
      titulo='{\delta}, X-Coherence ';
    case 2
      fifi='rawCohMnThNarrow';
      ppos=[.33+.33/2+xmarg .8/4*2+ymarg .33/2-2*xmarg .8/4-2*ymarg];
      titulo='{\theta}_{narrow}, X-Coherence ';
    case 3
      fifi='rawCohMnGa';
      ppos=[.33+.33/2+xmarg .8/4*1+ymarg .33/2-2*xmarg .8/4-2*ymarg];
      titulo='{\gamma}, X-Coherence ';
    case 4
      fifi='rawCohMnRi';
      ppos=[.33+.33/2+xmarg .8/4*0+ymarg .33/2-2*xmarg .8/4-2*ymarg];
      titulo='Rip, X-Coherence ';
  end
  ylab=' ';
  if isfield(r,fifi)
    subplot('position',ppos);
    hold on
    for i=bpix
      % assume that if theta ceherence exists it had also been calculated in other bands
      if length(diag(r(i).rawCohMnThNarrow))>1
        % more handy constants..
        pcol=AP.segmentType{i,3}; psym=AP.segmentType{i,4};
        if sigi==1
          tmpMn=cat(1,r(i).rawCohMnDe{1:LFPpcInd2,LFPpcInd2});
          tmpMn=[tmpMn; cat(1,r(i).rawCohMnDe{LFPpcInd2,LFPpcInd2+1:end})];
          tmpStd=cat(1,r(i).rawCohStdDe{1:LFPpcInd2,LFPpcInd2});
          tmpStd=[tmpStd; cat(1,r(i).rawCohStdDe{LFPpcInd2,LFPpcInd2+1:end})];
        elseif sigi==2
          tmpMn=cat(1,r(i).rawCohMnThNarrow{1:LFPpcInd2,LFPpcInd2});
          tmpMn=[tmpMn; cat(1,r(i).rawCohMnThNarrow{LFPpcInd2,LFPpcInd2+1:end})];
          tmpStd=cat(1,r(i).rawCohStdThNarrow{1:LFPpcInd2,LFPpcInd2});
          tmpStd=[tmpStd; cat(1,r(i).rawCohStdThNarrow{LFPpcInd2,LFPpcInd2+1:end})];
        elseif sigi==3
          tmpMn=cat(1,r(i).rawCohMnGa{1:LFPpcInd2,LFPpcInd2});
          tmpMn=[tmpMn; cat(1,r(i).rawCohMnGa{LFPpcInd2,LFPpcInd2+1:end})];
          tmpStd=cat(1,r(i).rawCohStdGa{1:LFPpcInd2,LFPpcInd2});
          tmpStd=[tmpStd; cat(1,r(i).rawCohStdGa{LFPpcInd2,LFPpcInd2+1:end})];
        else
          tmpMn=cat(1,r(i).rawCohMnRi{1:LFPpcInd2,LFPpcInd2});
          tmpMn=[tmpMn; cat(1,r(i).rawCohMnRi{LFPpcInd2,LFPpcInd2+1:end})];
          tmpStd=cat(1,r(i).rawCohStdRi{1:LFPpcInd2,LFPpcInd2});
          tmpStd=[tmpStd; cat(1,r(i).rawCohStdRi{LFPpcInd2,LFPpcInd2+1:end})];
        end
        % errorbars suspended until meaningful quantity available
        % ph=errorbar(WP.elx,tmpMn,zeros(size(tmpMn)),tmpStd,[psym '-']);
        ph=plot(WP.elx,tmpMn,[psym '-']);
        set(ph,'color',pcol);
      end
    end
    nicexax;
    set(gca,'ylim',[0 1.05],'xtick',WP.elx,'xticklabel',WP.xtl1);
    ylabel(ylab);
    % line to indicate principal electrode
    lh=line(WP.elx(LFPpcInd2)*[1 1]',get(gca,'ylim'));
    set(lh,'color',pcCol,'linestyle',':');
    % push to background
    set(gca,'children',circshift(get(gca,'children'),-1))
    % second axis for electrode #
    ax2=axes('position',get(gca,'position'),'color','none','xaxisloc','top',...
      'yaxisloc','right','ylim',get(gca,'ylim'),'ytick',get(gca,'ytick'),...
      'xlim',get(gca,'xlim'),'xtick',WP.elx,'xticklabel',WP.xtl2,'yticklabel',[]);
    th=title(titulo); set(th,'fontsize',10,'fontweight','bold');
    tp=get(th,'position');
    yl=get(gca,'ylim');
    tp(2)=yl(1)+diff(yl)*1.085;
    set(th,'position',tp);
  end
end


% ------ theta peak in power spec: frequency
if isfield(r,'rawPPeakTMn')
  subplot('position',[.66+xmarg .8/3*2+ymarg .33/2-2*xmarg .8/3-2*ymarg]);
  hold on
  for i=bpix
    if length(diag(r(i).rawPPeakTMn))>1
      % more handy constants..
      pcol=AP.segmentType{i,3}; psym=AP.segmentType{i,4};
      % nans will be ignored on plot
      ph=errorbar(WP.elx,cat(2,r(i).rawPPeakTMn{AP.dixie}),cat(2,r(i).rawPPeakTStd{AP.dixie}),[psym '-']);
      % ph=plot(WP.elx,cat(2,r(i).rawPPeakTMn{AP.dixie}),[psym '-']);
      set(ph,'color',pcol);
    end
  end
  nicexax;
  set(gca,'ylim',AP.peakFRng,'xtick',WP.elx,'xticklabel',WP.xtl1);
  ylabel('frequency (Hz)');
  % line to indicate principal electrode
  lh=line(WP.elx(LFPpcInd2)*[1 1]',get(gca,'ylim'));
  set(lh,'color',pcCol,'linestyle',':');
  % push to background
  set(gca,'children',circshift(get(gca,'children'),-1))
  % second axis for electrode #
  ax2=axes('position',get(gca,'position'),'color','none','xaxisloc','top',...
    'yaxisloc','right','ylim',get(gca,'ylim'),'ytick',get(gca,'ytick'),...
    'xlim',get(gca,'xlim'),'xtick',WP.elx,'xticklabel',WP.xtl2,'yticklabel',[]);
  th=title('{\theta}, Spec: F_{peak}'); set(th,'fontsize',10,'fontweight','bold');
  tp=get(th,'position');
  yl=get(gca,'ylim');
  tp(2)=yl(1)+diff(yl)*1.085;
  set(th,'position',tp);
end

% ------ theta peak in power spec: amplitude
if isfield(r,'rawPPeakMn')
  subplot('position',[.66+xmarg .8/3+ymarg .33/2-2*xmarg .8/3-2*ymarg]);
  hold on
  for i=bpix
    if length(diag(r(i).rawPPeakMn))>1
      % more handy constants..
      pcol=AP.segmentType{i,3}; psym=AP.segmentType{i,4};
      % nans will be ignored on plot
      % ph=plot(WP.elx,cat(2,r(i).rawPMnPeak{AP.dixie}),[psym '-']);
      ph=errorbar(WP.elx,cat(2,r(i).rawPPeakMn{AP.dixie}),cat(2,r(i).rawPPeakStd{AP.dixie}),[psym '-']);
      set(ph,'color',pcol);
    end
  end
  % line (below) will not be drawn if lower y limit is <0 and then yscale set to
  % log, so set y axis limits manually
  nicexax;
  yl=get(gca,'ylim');
  yl=yl.*[.8 1.2];
  set(gca,'ylim',yl);
  % set(gca,'yscale','log','xtick',WP.elx,'xticklabel',WP.xtl1);
  set(gca,'xtick',WP.elx,'xticklabel',WP.xtl1);
  ylabel('peak psd (mV^2/Hz)');
  % line to indicate principal electrode
  lh=line(WP.elx(LFPpcInd2)*[1 1]',get(gca,'ylim'));
  set(lh,'color',pcCol,'linestyle',':');
  % push to background
  set(gca,'children',circshift(get(gca,'children'),-1))
  % second axis for electrode #
  ax2=axes('position',get(gca,'position'),'color','none','xaxisloc','top',...
    'yaxisloc','right','ylim',get(gca,'ylim'),'ytick',get(gca,'ytick'),...
    'xlim',get(gca,'xlim'),'xtick',WP.elx,'xticklabel',WP.xtl2,'yticklabel',[]);
  th=title('{\theta}, Spec: A_{peak}'); set(th,'fontsize',10,'fontweight','bold');
  tp=get(th,'position');
  yl=get(gca,'ylim');
  tp(2)=yl(1)+diff(yl)*1.085;
  set(th,'position',tp);
end

% ------ theta peak in CROSS spec: frequency
if isfield(r,'rawPPeakTMn')
  subplot('position',[.66+.33/2+xmarg .8/3*2+ymarg .33/2-2*xmarg .8/3-2*ymarg]);
  hold on
  for i=bpix
    if length(diag(r(i).rawPPeakTMn))>1
      % more handy constants..
      pcol=AP.segmentType{i,3}; psym=AP.segmentType{i,4};
      % nans will be ignored on plot
      tmpMn=cat(1,r(i).rawPPeakTMn{1:LFPpcInd2,LFPpcInd2});
      tmpMn=[tmpMn; cat(1,r(i).rawPPeakTMn{LFPpcInd2,LFPpcInd2+1:end})];
      tmpStd=cat(1,r(i).rawPPeakTStd{1:LFPpcInd2,LFPpcInd2});
      tmpStd=[tmpStd; cat(1,r(i).rawPPeakTStd{LFPpcInd2,LFPpcInd2+1:end})];
      % ph=plot(WP.elx,tmpMn,[psym '-']);
      ph=errorbar(WP.elx,tmpMn,tmpStd,[psym '-']);
      set(ph,'color',pcol);
    end
  end
  nicexax;
  set(gca,'ylim',AP.peakFRng,'xtick',WP.elx,'xticklabel',WP.xtl1);
  ylabel('frequency (Hz)');
  % line to indicate principal electrode
  lh=line(WP.elx(LFPpcInd2)*[1 1]',get(gca,'ylim'));
  set(lh,'color',pcCol,'linestyle',':');
  % push to background
  set(gca,'children',circshift(get(gca,'children'),-1))
  % second axis for electrode #
  ax2=axes('position',get(gca,'position'),'color','none','xaxisloc','top',...
    'yaxisloc','right','ylim',get(gca,'ylim'),'ytick',get(gca,'ytick'),...
    'xlim',get(gca,'xlim'),'xtick',WP.elx,'xticklabel',WP.xtl2,'yticklabel',[]);
  th=title('{\theta}, X-Spec: F_{peak}'); set(th,'fontsize',10,'fontweight','bold');
  tp=get(th,'position');
  yl=get(gca,'ylim');
  tp(2)=yl(1)+diff(yl)*1.085;
  set(th,'position',tp);
end

% ------ theta peak in CROSS spec: amplitude
if isfield(r,'rawPPeakMn')
  subplot('position',[.66+.33/2+xmarg .8/3+ymarg .33/2-2*xmarg .8/3-2*ymarg]);
  hold on
  for i=bpix
    if length(diag(r(i).rawPPeakMn))>1
      % more handy constants..
      pcol=AP.segmentType{i,3}; psym=AP.segmentType{i,4};
      % nans will be ignored on plot
      tmpMn=cat(1,r(i).rawPPeakMn{1:LFPpcInd2,LFPpcInd2});
      tmpMn=[tmpMn; cat(1,r(i).rawPPeakMn{LFPpcInd2,LFPpcInd2+1:end})];
      tmpStd=cat(1,r(i).rawPPeakStd{1:LFPpcInd2,LFPpcInd2});
      tmpStd=[tmpStd; cat(1,r(i).rawPPeakStd{LFPpcInd2,LFPpcInd2+1:end})];
      % ph=plot(WP.elx,tmpMn,[psym '-']);
      ph=errorbar(WP.elx,tmpMn,tmpStd,[psym '-']);
      set(ph,'color',pcol);
    end
  end
  nicexyax;
  % set(gca,'yscale','log','xtick',WP.elx,'xticklabel',WP.xtl1);
  set(gca,'xtick',WP.elx,'xticklabel',WP.xtl1);
  ylabel('peak psd (mV^2/Hz)');
  % line to indicate principal electrode
  lh=line(WP.elx(LFPpcInd2)*[1 1]',get(gca,'ylim'));
  set(lh,'color',pcCol,'linestyle',':');
  % push to background
  set(gca,'children',circshift(get(gca,'children'),-1))
  % second axis for electrode #
  ax2=axes('position',get(gca,'position'),'color','none','xaxisloc','top',...
    'yaxisloc','right','ylim',get(gca,'ylim'),'ytick',get(gca,'ytick'),...
    'xlim',get(gca,'xlim'),'xtick',WP.elx,'xticklabel',WP.xtl2,'yticklabel',[]);
  th=title('{\theta}, X-Spec: A_{peak}'); set(th,'fontsize',10,'fontweight','bold');
  tp=get(th,'position');
  yl=get(gca,'ylim');
  tp(2)=yl(1)+diff(yl)*1.085;
  set(th,'position',tp);
end

% print this summary figure?
if ~isempty(AP.printas{1}),
  for i=1:length(AP.printas)
    pa=AP.printas{i};
    if strfind(pa,'ps'), ext='.ps';
    elseif strfind(pa,'jpeg'), ext='.jpg';
    else ext='';
    end
    print(pa,[figName '_specSum' ext]);
  end
end

% -------------------------------------------------------------------------------
%    ***  FIGURE IIa: summary plots of correlation analyses, first set ***
% -------------------------------------------------------------------------------
fh2a=mkfig('CCsum');
orient tall
labelscale('fontSz',7,'scaleFac',1.0,'lineW',1.0,'markSz',4);
% title
subplot('position',[xmarg .95 1-2*xmarg .02]);
set(gca,'xlim',[0 1],'ylim',[0 1]);
axis off
th=text(.5,.5,figTitle,'fontsize',12,'fontweight','bold','horizontalalignment','center');
drawnow
% # of rows = # of CC types
nRows=5;
% # of cols = # of result categories (currently: amplitude, lag, peak decay, Z-score(amp))
nCols=4;
xmarg=.04;
ymarg=.035;
xpos=(0:nCols-1)*((1-2*xmarg)/nCols)+2*xmarg;
xlen=(1-2*xmarg)/nCols-2*xmarg;
% leave additional space (.05)for title subplot
ypos=fliplr((0:nRows-1)*((1-.05-2*ymarg)/nRows)+2*ymarg);
ylen=(1-.05-2*ymarg)/nRows-2*ymarg;

% ------ delta CC:
row=1; lct=1;
if isfield(r,'deCCMn')
  % -- peaks
  sph=subplot('position',[xpos(lct) ypos(row) xlen ylen]);
  hold on
  for i=bpix
    if length(diag(r(i).deCCMn))>1
      % more handy constants..
      pcol=AP.segmentType{i,3}; psym=AP.segmentType{i,4};
      % extract CC data for principal channel: non analyzed channels
      % (nans) will be ignored on the plot
      tmpMn=cat(1, r(i).deCCPeakMn{1:LFPpcInd2,LFPpcInd2});
      tmpMn=[tmpMn; cat(1, r(i).deCCPeakMn{LFPpcInd2,LFPpcInd2+1:end})];
      tmpStd=cat(1, r(i).deCCPeakStd{1:LFPpcInd2,LFPpcInd2});
      tmpStd=[tmpStd; cat(1, r(i).deCCPeakStd{LFPpcInd2,LFPpcInd2+1:end})];
      ph=errorbar(WP.elx,tmpMn,tmpStd,[psym '-']);
      set(ph,'color',pcol);
    end
  end
  nicexax;
  set(gca,'ylim',[0 1.1],'xtick',WP.elx,'xticklabel',WP.xtl1);
  % line to indicate principal electrode
  lh=line(WP.elx(LFPpcInd2)*[1 1]',get(gca,'ylim'));
  set(lh,'color',pcCol,'linestyle',':');
  % push to background
  set(gca,'children',circshift(get(gca,'children'),-1))
  % inflate plot horizontally
  rexy('ax',gca,'xfac',1.2,'yfac',1.0);
  % second axis for electrode #
  ax2=axes('position',get(gca,'position'),'color','none','xaxisloc','top',...
    'yaxisloc','right','ylim',get(gca,'ylim'),'ytick',get(gca,'ytick'),...
    'xlim',get(gca,'xlim'),'xtick',WP.elx,'xticklabel',WP.xtl2,'yticklabel',[]);
  th=title('{\delta}, XCorr: A_{peak}'); set(th,'fontsize',10,'fontweight','bold');
  tp=get(th,'position');
  yl=get(gca,'ylim');
  tp(2)=yl(1)+diff(yl)*1.085;
  set(th,'position',tp);
  % -- lag
  lct=lct+1;
  sph=subplot('position',[xpos(lct) ypos(row) xlen ylen]);
  hold on
  for i=bpix
    if length(diag(r(i).deCCMn))>1
      pcol=AP.segmentType{i,3}; psym=AP.segmentType{i,4};
      tmpMn=cat(1, r(i).deCCPeakTMn{1:LFPpcInd2,LFPpcInd2});
      % symmetry - don't forget to invert
      tmpMn=[tmpMn; -1*cat(1, r(i).deCCPeakTMn{LFPpcInd2,LFPpcInd2+1:end})];
      tmpStd=cat(1, r(i).deCCPeakTStd{1:LFPpcInd2,LFPpcInd2});
      tmpStd=[tmpStd; cat(1, r(i).deCCPeakTStd{LFPpcInd2,LFPpcInd2+1:end})];
      ph=errorbar(WP.elx,tmpMn,tmpStd,[psym '-']);
      set(ph,'color',pcol);
    end
  end
  nicexyax;
  set(gca,'xtick',WP.elx,'xticklabel',WP.xtl1);
  ylabel('mean lag (ms)');
  % line to indicate principal electrode
  lh=line(WP.elx(LFPpcInd2)*[1 1]',get(gca,'ylim'));
  set(lh,'color',pcCol,'linestyle',':');
  % push to background
  set(gca,'children',circshift(get(gca,'children'),-1))
  % inflate plot horizontally
  rexy('ax',gca,'xfac',1.2,'yfac',1.0);
  % second axis for electrode #
  ax2=axes('position',get(gca,'position'),'color','none','xaxisloc','top',...
    'yaxisloc','right','ylim',get(gca,'ylim'),'ytick',get(gca,'ytick'),'xlim',get(gca,'xlim'),'xtick',WP.elx,'xticklabel',WP.xtl2,...
    'yticklabel',[]);
  th=title('{\delta}, XCorr: Lag_{peak}'); set(th,'fontsize',10,'fontweight','bold');
end

% ------ delta thetaEnv CC:
row=2; lct=1;
if isfield(r,'detheCCMn')
  % -- peaks
  sph=subplot('position',[xpos(lct) ypos(row) xlen ylen]);
  hold on
  for i=bpix
    if ~isempty(r(i).detheCCMn)
      % more handy constants..
      pcol=AP.segmentType{i,3}; psym=AP.segmentType{i,4};
      % for the sake of consistency with the plots above expand the arrays
      % to be plotted and fill the slots of non-computed channels with nans
      tmpMn=tempo;
      tmpStd=tempo;
      tmpMn(LFPccInd)=r(i).detheCCPeakMn;
      tmpStd(LFPccInd)=r(i).detheCCPeakStd;
      ph=errorbar(WP.elx,tmpMn,tmpStd,[psym '-']);
      set(ph,'color',pcol);
    end
  end
  nicexax;
  yl=get(gca,'ylim');
  if mean(yl)>0, yl=[0 .4];
  else yl=[-.4 0];
  end
  set(gca,'ylim',yl,'xtick',WP.elx,'xticklabel',WP.xtl1);
  % line to indicate principal electrode
  lh=line(WP.elx(LFPpcInd2)*[1 1]',get(gca,'ylim'));
  set(lh,'color',pcCol,'linestyle',':');
  % push to background
  set(gca,'children',circshift(get(gca,'children'),-1))
  % inflate plot horizontally
  rexy('ax',gca,'xfac',1.2,'yfac',1.0);
  % second axis for electrode #
  ax2=axes('position',get(gca,'position'),'color','none','xaxisloc','top','yaxisloc','right',...
    'ylim',get(gca,'ylim'),'ytick',get(gca,'ytick'),'xlim',get(gca,'xlim'),'xtick',WP.elx,...
    'xticklabel',WP.xtl2,'yticklabel',[]);
  th=title('{\delta}-{\theta}Env, Corr: A_{peak} '); set(th,'fontsize',10,'fontweight','bold');
  tp=get(th,'position');
  yl=get(gca,'ylim');
  tp(2)=yl(1)+diff(yl)*1.085;
  set(th,'position',tp);
  % -- lag
  lct=lct+1;
  sph=subplot('position',[xpos(lct) ypos(row) xlen ylen]);
  hold on
  for i=bpix
    if ~isempty(r(i).detheCCMn)
      pcol=AP.segmentType{i,3}; psym=AP.segmentType{i,4};
      % same as above applies
      tmpMn=tempo;
      tmpStd=tempo;
      tmpMn(LFPccInd)=r(i).detheCCPeakTMn;
      tmpStd(LFPccInd)=r(i).detheCCPeakTStd;
      ph=errorbar(WP.elx,tmpMn,tmpStd,[psym '-']);
      set(ph,'color',pcol);
    end
  end
  nicexyax;
  set(gca,'xtick',WP.elx,'xticklabel',WP.xtl1);
  ylabel('mean lag (ms)');
  % line to indicate principal electrode
  lh=line(WP.elx(LFPpcInd2)*[1 1]',get(gca,'ylim'));
  set(lh,'color',pcCol,'linestyle',':');
  % push to background
  set(gca,'children',circshift(get(gca,'children'),-1))
  % inflate plot horizontally
  rexy('ax',gca,'xfac',1.2,'yfac',1.0);
  % second axis for electrode #
  ax2=axes('position',get(gca,'position'),'color','none','xaxisloc','top','yaxisloc','right',...
    'ylim',get(gca,'ylim'),'ytick',get(gca,'ytick'),'xlim',get(gca,'xlim'),'xtick',WP.elx,...
    'xticklabel',WP.xtl2,'yticklabel',[]);
  th=title('{\delta}-{\theta}Env, Corr: Lag_{peak} '); set(th,'fontsize',10,'fontweight','bold');
end
% -- Z test stats
if isfield(r,'detheCCZTestP')
  lct=lct+1;
  sph=subplot('position',[xpos(lct) ypos(row) xlen ylen]);
  hold on
  for i=bpix
    if ~isempty(r(i).detheCCZTestP)
      pcol=AP.segmentType{i,3}; psym=AP.segmentType{i,4};
      % same as above applies
      tmpMn=tempo;
      tmpStd=tempo;
      tmpMn(LFPccInd)=mean(r(i).detheCCZTestP);
      tmpStd(LFPccInd)=std(r(i).detheCCZTestP);
      % curb absurdly significant values so we can see things in the interesting
      % range
      tmpMn(tmpMn<1e-6)=1e-6;
      tmpStd(tmpMn<1e-20)=nan;
      ph=errorbar(WP.elx,tmpMn,zeros(size(tmpMn)),tmpStd,[psym '-']);
      set(ph,'color',pcol);
    end
  end
  set(gca,'yscale','log');
  nicexax;
  yl=get(gca,'ylim');
  yl=[yl(1) 1];
  set(gca,'xtick',WP.elx,'xticklabel',WP.xtl1,'ylim',yl,'ytick',fliplr([10.^(0:-1:log10(yl(1)))]));
  % line to indicate the magical 0.05
  lh=line(get(gca,'xlim'),[.05 .05],'linestyle',':','color','r');
  % line to indicate principal electrode
  lh=line(WP.elx(LFPpcInd2)*[1 1]',get(gca,'ylim'));
  set(lh,'color',pcCol,'linestyle',':');
  % push to background
  set(gca,'children',circshift(get(gca,'children'),-1))
  % inflate plot horizontally
  rexy('ax',gca,'xfac',1.2,'yfac',1.0);
  % second axis for electrode #
  ax2=axes('position',get(gca,'position'),'color','none','xaxisloc','top','yaxisloc','right',...
    'ylim',get(gca,'ylim'),'ytick',get(gca,'ytick'),'xlim',get(gca,'xlim'),'xtick',WP.elx,...
    'xticklabel',WP.xtl2,'yticklabel',[]);
  th=title('{\delta}-{\theta}Env, Corr: p(Z-test)'); set(th,'fontsize',10,'fontweight','bold');
end

% ------ theta CC:
row=3; lct=1;
if isfield(r,'thCCMn')
  % -- peaks
  sph=subplot('position',[xpos(lct) ypos(row) xlen ylen]);
  lct=lct+1;
  hold on
  for i=bpix
    if length(diag(r(i).thCCMn))>1
      % more handy constants..
      pcol=AP.segmentType{i,3}; psym=AP.segmentType{i,4};
      % extract CC data for principal channel: non analyzed channels
      % (nans) will be ignored on the plot
      tmpMn=cat(1, r(i).thCCPeakMn{1:LFPpcInd2,LFPpcInd2});
      tmpMn=[tmpMn; cat(1, r(i).thCCPeakMn{LFPpcInd2,LFPpcInd2+1:end})];
      tmpStd=cat(1, r(i).thCCPeakStd{1:LFPpcInd2,LFPpcInd2});
      tmpStd=[tmpStd; cat(1, r(i).thCCPeakStd{LFPpcInd2,LFPpcInd2+1:end})];
      ph=errorbar(WP.elx,tmpMn,tmpStd,[psym '-']);
      set(ph,'color',pcol);
    end
  end
  nicexax;
  set(gca,'ylim',[0 1.1],'xtick',WP.elx,'xticklabel',WP.xtl1);
  % line to indicate principal electrode
  lh=line(WP.elx(LFPpcInd2)*[1 1]',get(gca,'ylim'));
  set(lh,'color',pcCol,'linestyle',':');
  % push to background
  set(gca,'children',circshift(get(gca,'children'),-1))
  % inflate plot horizontally
  rexy('ax',gca,'xfac',1.2,'yfac',1.0);
  % second axis for electrode #
  ax2=axes('position',get(gca,'position'),'color','none','xaxisloc','top',...
    'yaxisloc','right','ylim',get(gca,'ylim'),'ytick',get(gca,'ytick'),...
    'xlim',get(gca,'xlim'),'xtick',WP.elx,'xticklabel',WP.xtl2,'yticklabel',[]);
  th=title('{\theta}, XCorr: A_{peak}'); set(th,'fontsize',10,'fontweight','bold');
  tp=get(th,'position');
  yl=get(gca,'ylim');
  tp(2)=yl(1)+diff(yl)*1.085;
  set(th,'position',tp);
  % -- lag
  sph=subplot('position',[xpos(lct) ypos(row) xlen ylen]);
  hold on
  for i=bpix
    if length(diag(r(i).thCCMn))>1
      pcol=AP.segmentType{i,3}; psym=AP.segmentType{i,4};
      tmpMn=cat(1, r(i).thCCPeakTMn{1:LFPpcInd2,LFPpcInd2});
      % symmetry - don't forget to invert
      tmpMn=[tmpMn; -1*cat(1, r(i).thCCPeakTMn{LFPpcInd2,LFPpcInd2+1:end})];
      tmpStd=cat(1, r(i).thCCPeakTStd{1:LFPpcInd2,LFPpcInd2});
      tmpStd=[tmpStd; cat(1, r(i).thCCPeakTStd{LFPpcInd2,LFPpcInd2+1:end})];
      ph=errorbar(WP.elx,tmpMn,tmpStd,[psym '-']);
      set(ph,'color',pcol);
    end
  end
  nicexyax;
  set(gca,'xtick',WP.elx,'xticklabel',WP.xtl1);
  ylabel('mean lag (ms)');
  % line to indicate principal electrode
  lh=line(WP.elx(LFPpcInd2)*[1 1]',get(gca,'ylim'));
  set(lh,'color',pcCol,'linestyle',':');
  % push to background
  set(gca,'children',circshift(get(gca,'children'),-1))
  % inflate plot horizontally
  rexy('ax',gca,'xfac',1.2,'yfac',1.0);
  % second axis for electrode #
  ax2=axes('position',get(gca,'position'),'color','none','xaxisloc','top',...
    'yaxisloc','right','ylim',get(gca,'ylim'),'ytick',get(gca,'ytick'),'xlim',get(gca,'xlim'),'xtick',WP.elx,'xticklabel',WP.xtl2,...
    'yticklabel',[]);
  th=title('{\theta}, XCorr: Lag_{peak}'); set(th,'fontsize',10,'fontweight','bold');
  % -- CC peak decay
  lct=lct+1;
  sph=subplot('position',[xpos(lct) ypos(row) xlen ylen]);
  hold on
  for i=bpix
    if length(diag(r(i).thCCMn))>1
      pcol=AP.segmentType{i,3}; psym=AP.segmentType{i,4};
      tmpMn=cat(1, r(i).thCCPosPeakDecayMn{AP.dixie});
      tmpStd=cat(1, r(i).thCCPosPeakDecayStd{AP.dixie});
      ph=errorbar(WP.elx,tmpMn,tmpStd,[psym '-']);
      set(ph,'color',pcol);
    end
  end
  nicexy0ax;
  set(gca,'xtick',WP.elx,'xticklabel',WP.xtl1);
  ylabel('half-width (ms)');
  % line to indicate principal electrode
  lh=line(WP.elx(LFPpcInd2)*[1 1]',get(gca,'ylim'));
  set(lh,'color',pcCol,'linestyle',':');
  % push to background
  set(gca,'children',circshift(get(gca,'children'),-1))
  % inflate plot horizontally
  rexy('ax',gca,'xfac',1.2,'yfac',1.0);
  % second axis for electrode #
  ax2=axes('position',get(gca,'position'),'color','none','xaxisloc','top',...
    'yaxisloc','right','ylim',get(gca,'ylim'),'ytick',get(gca,'ytick'),'xlim',get(gca,'xlim'),'xtick',WP.elx,'xticklabel',WP.xtl2,...
    'yticklabel',[]);
  th=title('{\theta}, XCorr: AC decay_{peak}'); set(th,'fontsize',10,'fontweight','bold');
end

% ------ theta gammaEnv CC:
row=4; lct=1;
if isfield(r,'thgaeCCMn')
  % -- peaks
  sph=subplot('position',[xpos(lct) ypos(row) xlen ylen]);
  hold on
  for i=bpix
    if ~isempty(r(i).thgaeCCMn)
      % more handy constants..
      pcol=AP.segmentType{i,3}; psym=AP.segmentType{i,4};
      % for the sake of consistency with the plots above expand the arrays
      % to be plotted and fill the slots of non-computed channels with nans
      tmpMn=tempo;
      tmpStd=tempo;
      tmpMn(LFPccInd)=r(i).thgaeCCPeakMn;
      tmpStd(LFPccInd)=r(i).thgaeCCPeakStd;
      ph=errorbar(WP.elx,tmpMn,tmpStd,[psym '-']);
      set(ph,'color',pcol);
    end
  end
  nicexax;
  yl=get(gca,'ylim');
  if mean(yl)>0, yl=[0 .5];
  else yl=[-.5 0];
  end
  set(gca,'ylim',yl,'xtick',WP.elx,'xticklabel',WP.xtl1);
  % line to indicate principal electrode
  lh=line(WP.elx(LFPpcInd2)*[1 1]',get(gca,'ylim'));
  set(lh,'color',pcCol,'linestyle',':');
  % push to background
  set(gca,'children',circshift(get(gca,'children'),-1))
  % inflate plot horizontally
  rexy('ax',gca,'xfac',1.2,'yfac',1.0);
  % second axis for electrode #
  ax2=axes('position',get(gca,'position'),'color','none','xaxisloc','top','yaxisloc','right',...
    'ylim',get(gca,'ylim'),'ytick',get(gca,'ytick'),'xlim',get(gca,'xlim'),'xtick',WP.elx,...
    'xticklabel',WP.xtl2,'yticklabel',[]);
  th=title('{\theta}-{\gamma}Env, Corr: A_{peak} '); set(th,'fontsize',10,'fontweight','bold');
  tp=get(th,'position');
  yl=get(gca,'ylim');
  tp(2)=yl(1)+diff(yl)*1.085;
  set(th,'position',tp);
  % -- lag
  lct=lct+1;
  sph=subplot('position',[xpos(lct) ypos(row) xlen ylen]);
  hold on
  for i=bpix
    if ~isempty(r(i).thgaeCCMn)
      pcol=AP.segmentType{i,3}; psym=AP.segmentType{i,4};
      % same as above applies
      tmpMn=tempo;
      tmpStd=tempo;
      tmpMn(LFPccInd)=r(i).thgaeCCPeakTMn;
      tmpStd(LFPccInd)=r(i).thgaeCCPeakTStd;
      ph=errorbar(WP.elx,tmpMn,tmpStd,[psym '-']);
      set(ph,'color',pcol);
    end
  end
  nicexyax;
  set(gca,'xtick',WP.elx,'xticklabel',WP.xtl1);
  ylabel('mean lag (ms)');
  % line to indicate principal electrode
  lh=line(WP.elx(LFPpcInd2)*[1 1]',get(gca,'ylim'));
  set(lh,'color',pcCol,'linestyle',':');
  % push to background
  set(gca,'children',circshift(get(gca,'children'),-1))
  % inflate plot horizontally
  rexy('ax',gca,'xfac',1.2,'yfac',1.0);
  % second axis for electrode #
  ax2=axes('position',get(gca,'position'),'color','none','xaxisloc','top','yaxisloc','right',...
    'ylim',get(gca,'ylim'),'ytick',get(gca,'ytick'),'xlim',get(gca,'xlim'),'xtick',WP.elx,...
    'xticklabel',WP.xtl2,'yticklabel',[]);
  th=title('{\theta}-{\gamma}Env, Corr: Lag_{peak} '); set(th,'fontsize',10,'fontweight','bold');
  % -- CC peak decay
  lct=lct+1;
  sph=subplot('position',[xpos(lct) ypos(row) xlen ylen]);
  hold on
  for i=bpix
    if ~isempty(r(i).thgaeCCMn)
      % more handy constants..
      pcol=AP.segmentType{i,3}; psym=AP.segmentType{i,4};
      % for the sake of consistency with the plots above expand the arrays
      % to be plotted and fill the slots of non-computed channels with nans
      tmpMn=tempo;
      tmpStd=tempo;
      tmpMn(LFPccInd)=r(i).thgaeCCPosPeakDecayMn;
      tmpStd(LFPccInd)=r(i).thgaeCCPosPeakDecayStd;
      ph=errorbar(WP.elx,tmpMn,tmpStd,[psym '-']);
      set(ph,'color',pcol);
    end
  end
  nicexy0ax;
  set(gca,'xtick',WP.elx,'xticklabel',WP.xtl1);
  % line to indicate principal electrode
  lh=line(WP.elx(LFPpcInd2)*[1 1]',get(gca,'ylim'));
  set(lh,'color',pcCol,'linestyle',':');
  % push to background
  set(gca,'children',circshift(get(gca,'children'),-1))
  % inflate plot horizontally
  rexy('ax',gca,'xfac',1.2,'yfac',1.0);
  % second axis for electrode #
  ax2=axes('position',get(gca,'position'),'color','none','xaxisloc','top','yaxisloc','right',...
    'ylim',get(gca,'ylim'),'ytick',get(gca,'ytick'),'xlim',get(gca,'xlim'),'xtick',WP.elx,...
    'xticklabel',WP.xtl2,'yticklabel',[]);
  th=title('{\theta}-{\gamma}Env, AC: decay_{peak} '); set(th,'fontsize',10,'fontweight','bold');
  tp=get(th,'position');
  yl=get(gca,'ylim');
  tp(2)=yl(1)+diff(yl)*1.085;
  set(th,'position',tp);
end
% -- Z scores
if isfield(r,'thgaeCCZTestP')
  lct=lct+1;
  sph=subplot('position',[xpos(lct) ypos(row) xlen ylen]);
  hold on
  for i=bpix
    if ~isempty(r(i).thgaeCCZTestP)
      pcol=AP.segmentType{i,3}; psym=AP.segmentType{i,4};
      % same as above applies
      tmpMn=tempo;
      tmpStd=tempo;
      tmpMn(LFPccInd)=nanmean(r(i).thgaeCCZTestP);
      tmpStd(LFPccInd)=nanstd(r(i).thgaeCCZTestP);
      % curb absurdly significant values so we can see things in the interesting
      % range
      tmpMn(tmpMn<1e-6)=1e-6;
      tmpStd(tmpMn<1e-20)=nan;
      ph=errorbar(WP.elx,tmpMn,zeros(size(tmpMn)),tmpStd,[psym '-']);
      set(ph,'color',pcol);
    end
  end
  set(gca,'yscale','log');
  nicexax;
  yl=get(gca,'ylim');
  yl=[yl(1) 1];
  set(gca,'xtick',WP.elx,'xticklabel',WP.xtl1,'ylim',yl,'ytick',fliplr([10.^(0:-1:log10(yl(1)))]));
  % line to indicate the magical 0.05
  lh=line(get(gca,'xlim'),[.05 .05],'linestyle',':','color','r');
  % line to indicate principal electrode
  lh=line(WP.elx(LFPpcInd2)*[1 1]',get(gca,'ylim'));
  set(lh,'color',pcCol,'linestyle',':');
  % push to background
  set(gca,'children',circshift(get(gca,'children'),-1))
  % inflate plot horizontally
  rexy('ax',gca,'xfac',1.2,'yfac',1.0);
  % second axis for electrode #
  ax2=axes('position',get(gca,'position'),'color','none','xaxisloc','top','yaxisloc','right',...
    'ylim',get(gca,'ylim'),'ytick',get(gca,'ytick'),'xlim',get(gca,'xlim'),'xtick',WP.elx,...
    'xticklabel',WP.xtl2,'yticklabel',[]);
  th=title('{\theta}-{\gamma}Env, Corr: p(Z-test)'); set(th,'fontsize',10,'fontweight','bold');
end

% ------ gamma CC:
row=5; lct=1;
if isfield(r,'gaCCMn')
  % -- peaks
  sph=subplot('position',[xpos(lct) ypos(row) xlen ylen]);
  lct=lct+1;
  hold on
  for i=bpix
    if length(diag(r(i).gaCCMn))>1
      % more handy constants..
      pcol=AP.segmentType{i,3}; psym=AP.segmentType{i,4};
      % extract CC data for principal channel: non analyzed channels
      % (nans) will be ignored on the plot
      tmpMn=cat(1, r(i).gaCCPeakMn{1:LFPpcInd2,LFPpcInd2});
      tmpMn=[tmpMn; cat(1, r(i).gaCCPeakMn{LFPpcInd2,LFPpcInd2+1:end})];
      tmpStd=cat(1, r(i).gaCCPeakStd{1:LFPpcInd2,LFPpcInd2});
      tmpStd=[tmpStd; cat(1, r(i).gaCCPeakStd{LFPpcInd2,LFPpcInd2+1:end})];
      ph=errorbar(WP.elx,tmpMn,tmpStd,[psym '-']);
      set(ph,'color',pcol);
    end
  end
  nicexax;
  set(gca,'ylim',[0 1.1],'xtick',WP.elx,'xticklabel',WP.xtl1);
  % line to indicate principal electrode
  lh=line(WP.elx(LFPpcInd2)*[1 1]',get(gca,'ylim'));
  set(lh,'color',pcCol,'linestyle',':');
  % push to background
  set(gca,'children',circshift(get(gca,'children'),-1));
  % inflate plot horizontally
  rexy('ax',gca,'xfac',1.2,'yfac',1.0);
  % second axis for electrode #
  ax2=axes('position',get(gca,'position'),'color','none','xaxisloc','top','yaxisloc','right',...
    'ylim',get(gca,'ylim'),'ytick',get(gca,'ytick'),'xlim',get(gca,'xlim'),'xtick',WP.elx,...
    'xticklabel',WP.xtl2,'yticklabel',[]);
  th=title('{\gamma}, XCorr: A_{peak}'); set(th,'fontsize',10,'fontweight','bold');
  tp=get(th,'position');
  yl=get(gca,'ylim');
  tp(2)=yl(1)+diff(yl)*1.085;
  set(th,'position',tp);
  % -- lag
  sph=subplot('position',[xpos(lct) ypos(row) xlen ylen]);
  hold on
  for i=bpix
    if length(diag(r(i).gaCCMn))>1
      pcol=AP.segmentType{i,3}; psym=AP.segmentType{i,4};
      tmpMn=cat(1, r(i).gaCCPeakTMn{1:LFPpcInd2,LFPpcInd2});
      % symmetry - don't forget to invert
      tmpMn=[tmpMn; -1*cat(1, r(i).gaCCPeakTMn{LFPpcInd2,LFPpcInd2+1:end})];
      tmpStd=cat(1, r(i).gaCCPeakTStd{1:LFPpcInd2,LFPpcInd2});
      tmpStd=[tmpStd; cat(1, r(i).gaCCPeakTStd{LFPpcInd2,LFPpcInd2+1:end})];
      ph=errorbar(WP.elx,tmpMn,tmpStd,[psym '-']);
      set(ph,'color',pcol);
    end
  end
  nicexyax;
  set(gca,'xtick',WP.elx,'xticklabel',WP.xtl1);
  ylabel('mean lag (ms)');
  % line to indicate principal electrode
  lh=line(WP.elx(LFPpcInd2)*[1 1]',get(gca,'ylim'));
  set(lh,'color',pcCol,'linestyle',':');
  % push to background
  set(gca,'children',circshift(get(gca,'children'),-1))
  % inflate plot horizontally
  rexy('ax',gca,'xfac',1.2,'yfac',1.0);
  % second axis for electrode #
  ax2=axes('position',get(gca,'position'),'color','none','xaxisloc','top','yaxisloc','right',...
    'ylim',get(gca,'ylim'),'ytick',get(gca,'ytick'),'xlim',get(gca,'xlim'),'xtick',WP.elx,...
    'xticklabel',WP.xtl2,'yticklabel',[]);
  th=title('{\gamma}, XCorr: Lag_{peak}'); set(th,'fontsize',10,'fontweight','bold');
end

if length(get(gcf,'children'))<=1
  % this prevents empty figures from remaining on screen and getting printed
  delete(gcf);
else
  % print this summary figure?
  if ~isempty(AP.printas{1}),
    for i=1:length(AP.printas)
      pa=AP.printas{i};
      if strfind(pa,'ps'), ext='.ps';
      elseif strfind(pa,'jpeg'), ext='.jpg';
      else ext='';
      end
      print(pa,[figName '_ccSum' ext]);
    end
  end
end
% -------------------------------------------------------------------------------
%    ***  FIGURE IIb: summary plots of correlation analyses, second set ***
% -------------------------------------------------------------------------------
fh2b=mkfig('CCsum2');
orient tall
labelscale('fontSz',7,'scaleFac',1.0,'lineW',1.0,'markSz',4);
% title
subplot('position',[xmarg .95 1-2*xmarg .02]);
set(gca,'xlim',[0 1],'ylim',[0 1]);
axis off
th=text(.5,.5,figTitle,'fontsize',12,'fontweight','bold','horizontalalignment','center');
drawnow
% # of rows = # of CC types
nRows=4;
% # of cols = # of result categories (currently: amplitude, lag, Env CC peak half-width)
% +1 for similar appearance with CC set I
nCols=4;
xmarg=.04;
ymarg=.035;
xpos=(0:nCols-1)*((1-2*xmarg)/nCols)+2*xmarg;
xlen=(1-2*xmarg)/nCols-2*xmarg;
% leave additional space (.05)for title subplot
ypos=fliplr((0:nRows-1)*((1-.05-2*ymarg)/nRows)+2*ymarg);
ylen=(1-.05-2*ymarg)/nRows-2*ymarg;

% ------ gamma env CC:
row=1; lct=1;
if isfield(r,'gaeCCMn')
  % -- peaks
  sph=subplot('position',[xpos(lct) ypos(row) xlen ylen]);
  hold on
  for i=bpix
    if length(diag(r(i).gaeCCMn))>1
      % more handy constants..
      pcol=AP.segmentType{i,3}; psym=AP.segmentType{i,4};
      % extract CC data for principal channel: non analyzed channels
      % (nans) will be ignored on the plot
      tmpMn=cat(1, r(i).gaeCCPeakMn{1:LFPpcInd2,LFPpcInd2});
      tmpMn=[tmpMn; cat(1, r(i).gaeCCPeakMn{LFPpcInd2,LFPpcInd2+1:end})];
      tmpStd=cat(1, r(i).gaeCCPeakStd{1:LFPpcInd2,LFPpcInd2});
      tmpStd=[tmpStd; cat(1, r(i).gaeCCPeakStd{LFPpcInd2,LFPpcInd2+1:end})];
      ph=errorbar(WP.elx,tmpMn,tmpStd,[psym '-']);
      set(ph,'color',pcol);
    end
  end
  nicexax;
  set(gca,'ylim',[0 1.1],'xtick',WP.elx,'xticklabel',WP.xtl1);
  % line to indicate principal electrode
  lh=line(WP.elx(LFPpcInd2)*[1 1]',get(gca,'ylim'));
  set(lh,'color',pcCol,'linestyle',':');
  % push to background
  set(gca,'children',circshift(get(gca,'children'),-1))
  % inflate plot horizontally
  rexy('ax',gca,'xfac',1.2,'yfac',1.0);
  % second axis for electrode #
  ax2=axes('position',get(gca,'position'),'color','none','xaxisloc','top',...
    'yaxisloc','right','ylim',get(gca,'ylim'),'ytick',get(gca,'ytick'),...
    'xlim',get(gca,'xlim'),'xtick',WP.elx,'xticklabel',WP.xtl2,'yticklabel',[]);
  th=title('{\gamma}Env, XCorr: A_{peak}'); set(th,'fontsize',10,'fontweight','bold');
  tp=get(th,'position');
  yl=get(gca,'ylim');
  tp(2)=yl(1)+diff(yl)*1.085;
  set(th,'position',tp);
  % -- lag
  lct=lct+1;
  sph=subplot('position',[xpos(lct) ypos(row) xlen ylen]);
  hold on
  for i=bpix
    if length(diag(r(i).gaeCCMn))>1
      pcol=AP.segmentType{i,3}; psym=AP.segmentType{i,4};
      tmpMn=cat(1, r(i).gaeCCPeakTMn{1:LFPpcInd2,LFPpcInd2});
      % symmetry - don't forget to invert
      tmpMn=[tmpMn; -1*cat(1, r(i).gaeCCPeakTMn{LFPpcInd2,LFPpcInd2+1:end})];
      tmpStd=cat(1, r(i).gaeCCPeakTStd{1:LFPpcInd2,LFPpcInd2});
      tmpStd=[tmpStd; cat(1, r(i).gaeCCPeakTStd{LFPpcInd2,LFPpcInd2+1:end})];
      ph=errorbar(WP.elx,tmpMn,tmpStd,[psym '-']);
      set(ph,'color',pcol);
    end
  end
  nicexyax;
  set(gca,'xtick',WP.elx,'xticklabel',WP.xtl1);
  ylabel('mean lag (ms)');
  % line to indicate principal electrode
  lh=line(WP.elx(LFPpcInd2)*[1 1]',get(gca,'ylim'));
  set(lh,'color',pcCol,'linestyle',':');
  % push to background
  set(gca,'children',circshift(get(gca,'children'),-1))
  % inflate plot horizontally
  rexy('ax',gca,'xfac',1.2,'yfac',1.0);
  % second axis for electrode #
  ax2=axes('position',get(gca,'position'),'color','none','xaxisloc','top',...
    'yaxisloc','right','ylim',get(gca,'ylim'),'ytick',get(gca,'ytick'),'xlim',get(gca,'xlim'),'xtick',WP.elx,'xticklabel',WP.xtl2,...
    'yticklabel',[]);
  th=title('{\gamma}Env, XCorr: Lag_{peak}'); set(th,'fontsize',10,'fontweight','bold');
end

% ------ thetaLo env CC:
row=2; lct=1;
if isfield(r,'thLoeCCMn')
  % -- peaks
  sph=subplot('position',[xpos(lct) ypos(row) xlen ylen]);
  hold on
  for i=bpix
    if length(diag(r(i).thLoeCCMn))>1
      % more handy constants..
      pcol=AP.segmentType{i,3}; psym=AP.segmentType{i,4};
      % extract CC data for principal channel: non analyzed channels
      % (nans) will be ignored on the plot
      tmpMn=cat(1, r(i).thLoeCCPeakMn{1:LFPpcInd2,LFPpcInd2});
      tmpMn=[tmpMn; cat(1, r(i).thLoeCCPeakMn{LFPpcInd2,LFPpcInd2+1:end})];
      tmpStd=cat(1, r(i).thLoeCCPeakStd{1:LFPpcInd2,LFPpcInd2});
      tmpStd=[tmpStd; cat(1, r(i).thLoeCCPeakStd{LFPpcInd2,LFPpcInd2+1:end})];
      ph=errorbar(WP.elx,tmpMn,tmpStd,[psym '-']);
      set(ph,'color',pcol);
    end
  end
  nicexax;
  set(gca,'ylim',[0 1.1],'xtick',WP.elx,'xticklabel',WP.xtl1);
  % line to indicate principal electrode
  lh=line(WP.elx(LFPpcInd2)*[1 1]',get(gca,'ylim'));
  set(lh,'color',pcCol,'linestyle',':');
  % push to background
  set(gca,'children',circshift(get(gca,'children'),-1))
  % inflate plot horizontally
  rexy('ax',gca,'xfac',1.2,'yfac',1.0);
  % second axis for electrode #
  ax2=axes('position',get(gca,'position'),'color','none','xaxisloc','top',...
    'yaxisloc','right','ylim',get(gca,'ylim'),'ytick',get(gca,'ytick'),...
    'xlim',get(gca,'xlim'),'xtick',WP.elx,'xticklabel',WP.xtl2,'yticklabel',[]);
  th=title('{\theta}LoEnv, XCorr: A_{peak}'); set(th,'fontsize',10,'fontweight','bold');
  tp=get(th,'position');
  yl=get(gca,'ylim');
  tp(2)=yl(1)+diff(yl)*1.085;
  set(th,'position',tp);
  % -- lag
  lct=lct+1;
  sph=subplot('position',[xpos(lct) ypos(row) xlen ylen]);
  hold on
  for i=bpix
    if length(diag(r(i).thLoeCCMn))>1
      pcol=AP.segmentType{i,3}; psym=AP.segmentType{i,4};
      tmpMn=cat(1, r(i).thLoeCCPeakTMn{1:LFPpcInd2,LFPpcInd2});
      % symmetry - don't forget to invert
      tmpMn=[tmpMn; -1*cat(1, r(i).thLoeCCPeakTMn{LFPpcInd2,LFPpcInd2+1:end})];
      tmpStd=cat(1, r(i).thLoeCCPeakTStd{1:LFPpcInd2,LFPpcInd2});
      tmpStd=[tmpStd; cat(1, r(i).thLoeCCPeakTStd{LFPpcInd2,LFPpcInd2+1:end})];
      ph=errorbar(WP.elx,tmpMn,tmpStd,[psym '-']);
      set(ph,'color',pcol);
    end
  end
  nicexyax;
  set(gca,'xtick',WP.elx,'xticklabel',WP.xtl1);
  ylabel('mean lag (ms)');
  % line to indicate principal electrode
  lh=line(WP.elx(LFPpcInd2)*[1 1]',get(gca,'ylim'));
  set(lh,'color',pcCol,'linestyle',':');
  % push to background
  set(gca,'children',circshift(get(gca,'children'),-1))
  % inflate plot horizontally
  rexy('ax',gca,'xfac',1.2,'yfac',1.0);
  % second axis for electrode #
  ax2=axes('position',get(gca,'position'),'color','none','xaxisloc','top',...
    'yaxisloc','right','ylim',get(gca,'ylim'),'ytick',get(gca,'ytick'),'xlim',get(gca,'xlim'),'xtick',WP.elx,'xticklabel',WP.xtl2,...
    'yticklabel',[]);
  th=title('{\theta}LoEnv, XCorr: Lag_{peak}'); set(th,'fontsize',10,'fontweight','bold');
end

% ------ thetaHi env CC:
row=3; lct=1;
if isfield(r,'thHieCCMn')
  % -- peaks
  sph=subplot('position',[xpos(lct) ypos(row) xlen ylen]);
  hold on
  for i=bpix
    if length(diag(r(i).thHieCCMn))>1
      % more handy constants..
      pcol=AP.segmentType{i,3}; psym=AP.segmentType{i,4};
      % extract CC data for principal channel: non analyzed channels
      % (nans) will be ignored on the plot
      tmpMn=cat(1, r(i).thHieCCPeakMn{1:LFPpcInd2,LFPpcInd2});
      tmpMn=[tmpMn; cat(1, r(i).thHieCCPeakMn{LFPpcInd2,LFPpcInd2+1:end})];
      tmpStd=cat(1, r(i).thHieCCPeakStd{1:LFPpcInd2,LFPpcInd2});
      tmpStd=[tmpStd; cat(1, r(i).thHieCCPeakStd{LFPpcInd2,LFPpcInd2+1:end})];
      ph=errorbar(WP.elx,tmpMn,tmpStd,[psym '-']);
      set(ph,'color',pcol);
    end
  end
  nicexax;
  set(gca,'ylim',[0 1.1],'xtick',WP.elx,'xticklabel',WP.xtl1);
  % line to indicate principal electrode
  lh=line(WP.elx(LFPpcInd2)*[1 1]',get(gca,'ylim'));
  set(lh,'color',pcCol,'linestyle',':');
  % push to background
  set(gca,'children',circshift(get(gca,'children'),-1))
  % inflate plot horizontally
  rexy('ax',gca,'xfac',1.2,'yfac',1.0);
  % second axis for electrode #
  ax2=axes('position',get(gca,'position'),'color','none','xaxisloc','top',...
    'yaxisloc','right','ylim',get(gca,'ylim'),'ytick',get(gca,'ytick'),...
    'xlim',get(gca,'xlim'),'xtick',WP.elx,'xticklabel',WP.xtl2,'yticklabel',[]);
  th=title('{\theta}HiEnv, XCorr: A_{peak}'); set(th,'fontsize',10,'fontweight','bold');
  tp=get(th,'position');
  yl=get(gca,'ylim');
  tp(2)=yl(1)+diff(yl)*1.085;
  set(th,'position',tp);
  % -- lag
  lct=lct+1;
  sph=subplot('position',[xpos(lct) ypos(row) xlen ylen]);
  hold on
  for i=bpix
    if length(diag(r(i).thHieCCMn))>1
      pcol=AP.segmentType{i,3}; psym=AP.segmentType{i,4};
      tmpMn=cat(1, r(i).thHieCCPeakTMn{1:LFPpcInd2,LFPpcInd2});
      % symmetry - don't forget to invert
      tmpMn=[tmpMn; -1*cat(1, r(i).thHieCCPeakTMn{LFPpcInd2,LFPpcInd2+1:end})];
      tmpStd=cat(1, r(i).thHieCCPeakTStd{1:LFPpcInd2,LFPpcInd2});
      tmpStd=[tmpStd; cat(1, r(i).thHieCCPeakTStd{LFPpcInd2,LFPpcInd2+1:end})];
      ph=errorbar(WP.elx,tmpMn,tmpStd,[psym '-']);
      set(ph,'color',pcol);
    end
  end
  nicexyax;
  set(gca,'xtick',WP.elx,'xticklabel',WP.xtl1);
  ylabel('mean lag (ms)');
  % line to indicate principal electrode
  lh=line(WP.elx(LFPpcInd2)*[1 1]',get(gca,'ylim'));
  set(lh,'color',pcCol,'linestyle',':');
  % push to background
  set(gca,'children',circshift(get(gca,'children'),-1))
  % inflate plot horizontally
  rexy('ax',gca,'xfac',1.2,'yfac',1.0);
  % second axis for electrode #
  ax2=axes('position',get(gca,'position'),'color','none','xaxisloc','top',...
    'yaxisloc','right','ylim',get(gca,'ylim'),'ytick',get(gca,'ytick'),'xlim',get(gca,'xlim'),'xtick',WP.elx,'xticklabel',WP.xtl2,...
    'yticklabel',[]);
  th=title('{\theta}HiEnv, XCorr: Lag_{peak}'); set(th,'fontsize',10,'fontweight','bold');
end

if length(get(gcf,'children'))<=1
  % this prevents empty figures from remaining on screen and getting printed
  delete(gcf);
else
  % print this summary figure?
  if ~isempty(AP.printas{1}),
    for i=1:length(AP.printas)
      pa=AP.printas{i};
      if strfind(pa,'ps'), ext='.ps';
      elseif strfind(pa,'jpeg'), ext='.jpg';
      else ext='';
      end
      print(pa,[figName '_ccSum2' ext]);
    end
  end
end

% -------------------------------------------------------------------------------
%               ***  FIGURE III: contour plots from spectral analysis ***
% -------------------------------------------------------------------------------
if spec_contourP
  fh3=mkfig('CohContour');
  orient tall
  labelscale('fontSz',10,'scaleFac',1.0,'lineW',1.0,'markSz',4);
  colormap(hot);
  % # of contours
  nLev=20;
  % # of rows = currently use 2 of 3 (for different F axes)
  nRows=length(bpix);
  % # of cols = # of behavior types with results
  nCols=2;
  xmarg=.04;
  ymarg=.035;
  xpos=(0:nCols-1)*((1-2*xmarg)/nCols)+2*xmarg;
  xlen=(1-2*xmarg)/nCols-2*xmarg;
  % leave additional space (.05)for title subplot
  ypos=fliplr((0:nRows-1)*((1-.05-2*ymarg)/nRows)+2*ymarg);
  ylen=(1-.05-2*ymarg)/nRows-2*ymarg;
  % title
  subplot('position',[xmarg .95 1-2*xmarg .02]);
  set(gca,'xlim',[0 1],'ylim',[0 1]);
  axis off
  th=text(.5,.5,figTitle,'fontsize',12,'fontweight','bold','horizontalalignment','center');
  if isfield(r,'rawCohMn')
    srfc=[];
    % large range (0-100 Hz)
    frix=find(r(1).F>=0 & r(1).F<=100);
    % around theta
    sfrix=find(r(1).F(frix)<=18);
    for fri=1:2
      % I.+II. coherence principal vs others (2 frequency ranges)
      row=fri; lct=1;
      hold on
      for i=bpix
        if length(diag(r(i).rawCohMn))>1
          % freq range from which average NARROW theta coherence had been calculated
          narrThF=r(i).rawPMnPeakT{LFPpcInd2,LFPpcInd2}+[-1 1];
          sph=subplot('position',[xpos(lct) ypos(row) xlen ylen]);
          % first run: generate surface array
          if fri==1
            srfc=cat(3, srfc, repmat(nan,length(frix),nAllLFPCh));
            % the lines below take into account non-analyzed channels (which have
            % single nans in the corresponding cells)
            for g=1:nAllLFPCh
              ix1=min(g,LFPpcInd2);
              ix2=max(g,LFPpcInd2);
              if length(r(i).rawCohMn{ix1,ix2})>1
                srfc(1:length(frix),g,lct)=r(i).rawCohMn{ix1,ix2}(frix);
              else
                srfc(1:length(frix),g,lct)=nan;
              end
            end
            % 2nd run: cut down
          elseif fri==2
            if lct==1
              % cut down original index and use it henceforth
              frix=frix(sfrix);
              srfc=srfc(frix,:,:);
            end
          end
          % --- contour plot
          [c,cph]=contourf(1:size(srfc,2),r(1).F(frix),srfc(:,:,lct),nLev);
          clear c;
          axis tight
          set(cph(:),'linestyle','none');
          cph=gca;
          % matlab V. 7.0.1 completely fucks up overlaid axes with colorbars, so let it be
          % for now
          % cbh=colorbar(cbPos);
          % use the contour plot's x axis as electrode # indicator
          set(cph,'xaxisloc','top','yaxisloc','right','ylim',[r(1).F(frix(1)) r(1).F(frix(end))],'ytick',[],...
            'xlim',[1 nAllLFPCh],'xtick',[1:nAllLFPCh],'xticklabel',WP.xtl2);
          % x axis limits must be fixed for aligned overlay with surface plot
          ax2=axes('position',get(cph,'position'),'color','none');
          hold on;
          set(gca,'ylim',[r(1).F(frix(1)) r(1).F(frix(end))],'xlim',[WP.elx(1) WP.elx(end)],...
            'xtick',WP.elx,'xticklabel',WP.xtl1);
          axis manual;
          pcol=AP.segmentType{i,3}; psym=AP.segmentType{i,4};
          grid on;
          % line to indicate principal electrode
          lh=line(WP.elx(LFPpcInd2)*[1 1]',get(gca,'ylim'));
          set(lh,'color',pcCol,'linestyle',':');
          % push to background
          set(gca,'children',circshift(get(gca,'children'),-1))
          % lines to indicate freq ranges from which average coherence had been calculated
          lh=line([WP.elx(1) WP.elx(end)]'*[1 1],[1;1]*AP.delta,'color','c');
          lh=line([WP.elx(1) WP.elx(end)]'*[1 1],[1;1]*narrThF,'color','b');
          lh=line([WP.elx(1) WP.elx(end)]'*[1 1],[1;1]*AP.gamma,'color','g');
          title(['Coherence;' r(i).segmentType]);
          if lct==1, ylabel('Frequency (Hz)'); end
          grid on;
          lct=lct+1;
        end
      end
    end
  end

  if length(get(gcf,'children'))<=1
    % this prevents empty figures from remaining on screen and getting printed
    delete(gcf);
  else
    % print this summary figure?
    if ~isempty(AP.printas{1}),
      for i=1:length(AP.printas)
        pa=AP.printas{i};
        if strfind(pa,'ps'), ext='.ps';
        elseif strfind(pa,'jpeg'), ext='.jpg';
        else ext='';
        end
        print(pa,[figName '_specContour' ext]);
      end
    end
  end
end
% -------------------------------------------------------------------------------
%    ***  FIGURE IVa: contour plots from correlation analyses, first set ***
% -------------------------------------------------------------------------------
if cc_contourP
  fh4a=mkfig('CCavContour');
  orient tall
  labelscale('fontSz',10,'scaleFac',1.0,'lineW',1.0,'markSz',4);
  colormap(bone);
  % # of contours
  nLev=20;
  % # of rows = # of CC types
  nRows=5;
  % # of cols = # of behavior types with results
  nCols=length(bpix);
  xmarg=.04;
  ymarg=.035;
  xpos=(0:nCols-1)*((1-2*xmarg)/nCols)+2*xmarg;
  xlen=(1-2*xmarg)/nCols-2*xmarg;
  % leave additional space (.05)for title subplot
  ypos=fliplr((0:nRows-1)*((1-.05-2*ymarg)/nRows)+2*ymarg);
  ylen=(1-.05-2*ymarg)/nRows-2*ymarg;
  % title
  subplot('position',[xmarg .95 1-2*xmarg .02]);
  set(gca,'xlim',[0 1],'ylim',[0 1]);
  axis off
  th=text(.5,.5,figTitle,'fontsize',12,'fontweight','bold','horizontalalignment','center');

  if isfield(r,'deCCMn')
    % I. delta (range to plot: the lesser of [+/- 500 ms, AP.ccLagPts])
    row=1; lct=1;
    ccw=[-1 1]*min(AP.ccLagPts,cont2discrete(500,osi*.001,'intv',1));
    cci=AP.ccLagPts+[ccw(1):ccw(2)]+1;
    ccp=discrete2cont(cci-AP.ccLagPts,osi*.001,'intv',0);
    hold on
    for i=bpix
      if length(diag(r(i).deCCMn))>1
        sph=subplot('position',[xpos(lct) ypos(row) xlen ylen]);
        lct=lct+1;
        srfc=repmat(nan,length(cci),nAllLFPCh);
        % the lines below take into account non-analyzed channels (which have
        % single nans in the corresponding cells)
        for g=1:nAllLFPCh
          ix1=min(g,LFPpcInd2);
          ix2=max(g,LFPpcInd2);
          if length(r(i).deCCMn{ix1,ix2})>1
            srfc(1:length(cci),g)=r(i).deCCMn{ix1,ix2}(cci);
            % symmetry
            if g>LFPpcInd2, srfc(:,g)=flipud(srfc(:,g)); end
          else
            srfc(1:length(cci),g)=nan;
          end
        end
        % --- contour plot
        [c,cph]=contourf(1:nAllLFPCh,ccp,srfc,nLev);
        clear c;
        axis tight
        set(cph(:),'linestyle','none');
        cph=gca;
        % cbh=colorbar(cbPos);
        % use the contour plot's x axis as electrode # indicator
        set(cph,'xaxisloc','top','yaxisloc','right','ytick',[],...
          'xlim',[1 nAllLFPCh],'xtick',[1:nAllLFPCh],'xticklabel',WP.xtl2);
        % x axis limits must be fixed for aligned overlay with surface plot
        ax2=axes('position',get(cph,'position'),'color','none');
        hold on;
        set(gca,'ylim',[ccp(1) ccp(end)],'xlim',[WP.elx(1) WP.elx(end)],...
          'xtick',WP.elx,'xticklabel',WP.xtl1);
        pcol=AP.segmentType{i,3}; psym=AP.segmentType{i,4};
        tmpMn=cat(1,r(i).deCCPeakTMn{1:LFPpcInd2,LFPpcInd2});
        % symmetry - don't forget to invert
        tmpMn=[tmpMn; -1*cat(1, r(i).deCCPeakTMn{LFPpcInd2,LFPpcInd2+1:end})];
        tmpStd=cat(1, r(i).deCCPeakTStd{1:LFPpcInd2,LFPpcInd2});
        tmpStd=[tmpStd; cat(1, r(i).deCCPeakTStd{LFPpcInd2,LFPpcInd2+1:end})];
        ph=errorbar(WP.elx,tmpMn,tmpStd,[psym '-']);
        set(ph,'color',pcol);
        if lct==2, ylabel('lag (ms)'); end
        grid on;
        % line to indicate principal electrode
        lh=line(WP.elx(LFPpcInd2)*[1 1]',get(gca,'ylim'));
        set(lh,'color',pcCol,'linestyle',':');
        % push to background
        set(gca,'children',circshift(get(gca,'children'),-1))
        title(['XC{\delta},{\delta}(ref);' r(i).segmentType]);
      end
    end
  end

  if isfield(r,'detheCCMn')
    % II. delta-thetaEnv (range to plot: the lesser of [+/- 500 ms, AP.ccLagPts])
    row=2; lct=1;
    ccw=[-1 1]*min(AP.ccLagPts,cont2discrete(500,osi*.001,'intv',1));
    cci=AP.ccLagPts+[ccw(1):ccw(2)]+1;
    ccp=discrete2cont(cci-AP.ccLagPts,osi*.001,'intv',0);
    % limits [-.33 .27]
    v=((0:nLev-1)/(nLev-1)-.55)*.1*6;
    hold on
    for i=bpix
      if ~isempty(r(i).detheCCMn)
        sph=subplot('position',[xpos(lct) ypos(row) xlen ylen]);
        lct=lct+1;
        srfc=repmat(nan,length(cci),nAllLFPCh);
        srfc(:,LFPccInd)=r(i).detheCCMn(cci,:);
        % --- contour plot
        [c,cph]=contourf(srfc,v);
        clear c;
        axis tight
        set(cph(:),'linestyle','none');
        set(gca,'clim',v([1 end]));
        cph=gca;
        % cbh=colorbar(cbPos);;
        % use the contour plot's x axis as electrode # indicator
        set(cph,'xaxisloc','top','yaxisloc','right','ytick',[],...
          'xlim',[1 nAllLFPCh],'xtick',[1:nAllLFPCh],'xticklabel',WP.xtl2);
        % x axis limits must be fixed for aligned overlay with surface plot
        ax2=axes('position',get(cph,'position'),'color','none');
        hold on;
        set(gca,'ylim',[ccp(1) ccp(end)],'xlim',[WP.elx(1) WP.elx(end)],...
          'xtick',WP.elx,'xticklabel',WP.xtl1);
        axis manual;
        pcol=AP.segmentType{i,3}; psym=AP.segmentType{i,4};
        tmpMn=tempo;
        tmpStd=tempo;
        tmpMn(LFPccInd)=r(i).detheCCPeakTMn;
        tmpStd(LFPccInd)=r(i).detheCCPeakTStd;
        ph=errorbar(WP.elx,tmpMn,tmpStd,[psym '-']);
        set(ph,'color',pcol);
        if lct==2, ylabel('lag (ms)'); end
        grid on;
        % line to indicate principal electrode
        lh=line(WP.elx(LFPpcInd2)*[1 1]',get(gca,'ylim'));
        set(lh,'color',pcCol,'linestyle',':');
        % push to background
        set(gca,'children',circshift(get(gca,'children'),-1))
        title(['XC{\delta},{\theta}Env;' r(i).segmentType]);
      end
    end
  end

  if isfield(r,'thCCMn')
    % III. theta (range to plot: the lesser of [+/- 200 ms, AP.ccLagPts])
    row=3; lct=1;
    ccw=[-1 1]*min(AP.ccLagPts,cont2discrete(200,osi*.001,'intv',1));
    cci=AP.ccLagPts+[ccw(1):ccw(2)]+1;
    ccp=discrete2cont(cci-AP.ccLagPts,osi*.001,'intv',0);
    hold on
    for i=bpix
      if length(diag(r(i).thCCMn))>1
        sph=subplot('position',[xpos(lct) ypos(row) xlen ylen]);
        lct=lct+1;
        srfc=repmat(nan,length(cci),nAllLFPCh);
        % the lines below take into account non-analyzed channels (which have
        % single nans in the corresponding cells)
        for g=1:nAllLFPCh
          ix1=min(g,LFPpcInd2);
          ix2=max(g,LFPpcInd2);
          if length(r(i).thCCMn{ix1,ix2})>1
            srfc(1:length(cci),g)=r(i).thCCMn{ix1,ix2}(cci);
            % symmetry
            if g>LFPpcInd2, srfc(:,g)=flipud(srfc(:,g)); end
          else
            srfc(1:length(cci),g)=nan;
          end
        end
        % --- contour plot
        [c,cph]=contourf(1:nAllLFPCh,ccp,srfc,nLev);
        clear c;
        axis tight
        set(cph(:),'linestyle','none');
        cph=gca;
        % cbh=colorbar(cbPos);
        % use the contour plot's x axis as electrode # indicator
        set(cph,'xaxisloc','top','yaxisloc','right','ytick',[],...
          'xlim',[1 nAllLFPCh],'xtick',[1:nAllLFPCh],'xticklabel',WP.xtl2);
        % x axis limits must be fixed for aligned overlay with surface plot
        ax2=axes('position',get(cph,'position'),'color','none');
        hold on;
        set(gca,'ylim',[ccp(1) ccp(end)],'xlim',[WP.elx(1) WP.elx(end)],...
          'xtick',WP.elx,'xticklabel',WP.xtl1);
        axis manual;
        pcol=AP.segmentType{i,3}; psym=AP.segmentType{i,4};
        tmpMn=cat(1,r(i).thCCPeakTMn{1:LFPpcInd2,LFPpcInd2});
        % symmetry - don't forget to invert
        tmpMn=[tmpMn; -1*cat(1, r(i).thCCPeakTMn{LFPpcInd2,LFPpcInd2+1:end})];
        tmpStd=cat(1, r(i).thCCPeakTStd{1:LFPpcInd2,LFPpcInd2});
        tmpStd=[tmpStd; cat(1, r(i).thCCPeakTStd{LFPpcInd2,LFPpcInd2+1:end})];
        ph=errorbar(WP.elx,tmpMn,tmpStd,[psym '-']);
        set(ph,'color',pcol);
        if lct==2, ylabel('lag (ms)'); end
        grid on;
        % line to indicate principal electrode
        lh=line(WP.elx(LFPpcInd2)*[1 1]',get(gca,'ylim'));
        set(lh,'color',pcCol,'linestyle',':');
        % push to background
        set(gca,'children',circshift(get(gca,'children'),-1))
        title(['XC{\theta},{\theta}(ref);' r(i).segmentType]);
      end
    end
  end

  if isfield(r,'thgaeCCMn')
    % IV. theta-gammaEnv (range to plot: the lesser of [+/- 1 slow theta period, AP.ccLagPts])
    row=4; lct=1;
    ccw=[-1 1]*min(AP.ccLagPts,cont2discrete(200,osi*.001,'intv',1));
    cci=AP.ccLagPts+[ccw(1):ccw(2)]+1;
    ccp=discrete2cont(cci-AP.ccLagPts,osi*.001,'intv',0);
    hold on
    for i=bpix
      if ~isempty(r(i).thgaeCCMn)
        sph=subplot('position',[xpos(lct) ypos(row) xlen ylen]);
        lct=lct+1;
        srfc=repmat(nan,length(cci),nAllLFPCh);
        srfc(:,LFPccInd)=r(i).thgaeCCMn(cci,:);
        % --- contour plot
        [c,cph]=contourf(1:nAllLFPCh,ccp,srfc,nLev);
        clear c;
        axis tight
        set(cph(:),'linestyle','none');
        cph=gca;
        % cbh=colorbar(cbPos);
        % use the contour plot's x axis as electrode # indicator
        set(cph,'xaxisloc','top','yaxisloc','right','ytick',[],...
          'xlim',[1 nAllLFPCh],'xtick',[1:nAllLFPCh],'xticklabel',WP.xtl2);
        % x axis limits must be fixed for aligned overlay with surface plot
        ax2=axes('position',get(cph,'position'),'color','none');
        hold on;
        set(gca,'ylim',[ccp(1) ccp(end)],'xlim',[WP.elx(1) WP.elx(end)],...
          'xtick',WP.elx,'xticklabel',WP.xtl1);
        axis manual;
        pcol=AP.segmentType{i,3}; psym=AP.segmentType{i,4};
        tmpMn=tempo;
        tmpStd=tempo;
        tmpMn(LFPccInd)=r(i).thgaeCCPeakTMn;
        tmpStd(LFPccInd)=r(i).thgaeCCPeakTStd;
        ph=errorbar(WP.elx,tmpMn,tmpStd,[psym '-']);
        set(ph,'color',pcol);
        if lct==2, ylabel('lag (ms)'); end
        grid on;
        % line to indicate principal electrode
        lh=line(WP.elx(LFPpcInd2)*[1 1]',get(gca,'ylim'));
        set(lh,'color',pcCol,'linestyle',':');
        % push to background
        set(gca,'children',circshift(get(gca,'children'),-1))
        title(['XC{\theta},{\gamma}Env;' r(i).segmentType]);
      end
    end
  end

  if isfield(r,'gaCCMn')
    % V. gamma (range to plot: the lesser of [+/- 1 slow gamma periods, AP.ccLagPts])
    row=5; lct=1;
    ccw=[-1 1]*min(AP.ccLagPts,cont2discrete(1*1000/AP.gamma(1),osi*.001,'intv',1));
    cci=AP.ccLagPts+[ccw(1):ccw(2)]+1;
    ccp=discrete2cont(cci-AP.ccLagPts,osi*.001,'intv',0);
    hold on
    for i=bpix
      if length(diag(r(i).gaCCMn))>1
        sph=subplot('position',[xpos(lct) ypos(row) xlen ylen]);
        lct=lct+1;
        srfc=repmat(nan,length(cci),nAllLFPCh);
        % the lines below take into account non-analyzed channels (which have
        % single nans in the corresponding cells)
        for g=1:nAllLFPCh
          ix1=min(g,LFPpcInd2);
          ix2=max(g,LFPpcInd2);
          if length(r(i).gaCCMn{ix1,ix2})>1
            srfc(1:length(cci),g)=r(i).gaCCMn{ix1,ix2}(cci);
            % symmetry
            if g>pcIdx, srfc(:,g)=flipud(srfc(:,g)); end
          else
            srfc(1:length(cci),g)=nan;
          end
        end
        % --- contour plot
        [c,cph]=contourf(1:nAllLFPCh,ccp,srfc,nLev);
        clear c;
        axis tight
        set(cph(:),'linestyle','none');
        cph=gca;
        % cbh=colorbar(cbPos);
        % use the contour plot's x axis as electrode # indicator
        set(cph,'xaxisloc','top','yaxisloc','right','ytick',[],...
          'xlim',[1 nAllLFPCh],'xtick',[1:nAllLFPCh],'xticklabel',WP.xtl2);
        % x axis limits must be fixed for aligned overlay with surface plot
        ax2=axes('position',get(cph,'position'),'color','none');
        hold on;
        set(gca,'ylim',[ccp(1) ccp(end)],'xlim',[WP.elx(1) WP.elx(end)],...
          'xtick',WP.elx,'xticklabel',WP.xtl1);
        axis manual;
        pcol=AP.segmentType{i,3}; psym=AP.segmentType{i,4};
        tmpMn=cat(1,r(i).gaCCPeakTMn{1:LFPpcInd2,LFPpcInd2});
        % symmetry - don't forget to invert
        tmpMn=[tmpMn; -1*cat(1, r(i).gaCCPeakTMn{LFPpcInd2,LFPpcInd2+1:end})];
        tmpStd=cat(1, r(i).gaCCPeakTStd{1:LFPpcInd2,LFPpcInd2});
        tmpStd=[tmpStd; cat(1, r(i).gaCCPeakTStd{LFPpcInd2,LFPpcInd2+1:end})];
        ph=errorbar(WP.elx,tmpMn,tmpStd,[psym '-']);
        set(ph,'color',pcol);
        if lct==2, ylabel('lag (ms)'); end
        grid on;
        % line to indicate principal electrode
        lh=line(WP.elx(LFPpcInd2)*[1 1]',get(gca,'ylim'));
        set(lh,'color',pcCol,'linestyle',':');
        % push to background
        set(gca,'children',circshift(get(gca,'children'),-1))
        title(['XC{\gamma},{\gamma}(ref);' r(i).segmentType]);
      end
    end
  end

  if length(get(gcf,'children'))<=1
    % this prevents empty figures from remaining on screen and getting printed
    delete(gcf);
  else
    % print this summary figure?
    if ~isempty(AP.printas{1}),
      for i=1:length(AP.printas)
        pa=AP.printas{i};
        if strfind(pa,'ps'), ext='.ps';
        elseif strfind(pa,'jpeg'), ext='.jpg';
        else ext='';
        end
        print(pa,[figName '_ccContour' ext]);
      end
    end
  end

  % -------------------------------------------------------------------------------
  %    ***  FIGURE IVb: contour plots from correlation analyses, second set ***
  % -------------------------------------------------------------------------------

  fh4b=mkfig('CCavContour2');
  orient tall
  labelscale('fontSz',10,'scaleFac',1.0,'lineW',1.0,'markSz',4);
  colormap(bone);
  % # of contours
  nLev=20;
  % # of rows = # of CC types
  nRows=4;
  % # of cols = # of behavior types with results
  nCols=length(bpix);
  xmarg=.04;
  ymarg=.035;
  xpos=(0:nCols-1)*((1-2*xmarg)/nCols)+2*xmarg;
  xlen=(1-2*xmarg)/nCols-2*xmarg;
  % leave additional space (.05)for title subplot
  ypos=fliplr((0:nRows-1)*((1-.05-2*ymarg)/nRows)+2*ymarg);
  ylen=(1-.05-2*ymarg)/nRows-2*ymarg;
  % title
  subplot('position',[xmarg .95 1-2*xmarg .02]);
  set(gca,'xlim',[0 1],'ylim',[0 1]);
  axis off
  th=text(.5,.5,figTitle,'fontsize',12,'fontweight','bold','horizontalalignment','center');

  if isfield(r,'gaeCCMn')
    % 1. gamma env (range to plot: the lesser of [+/- 2 slow gamma periods, AP.ccLagPts])
    row=1; lct=1;
    ccw=[-1 1]*min(AP.ccLagPts,cont2discrete(2*1000/AP.gamma(1),osi*.001,'intv',1));
    cci=AP.ccLagPts+[ccw(1):ccw(2)]+1;
    ccp=discrete2cont(cci-AP.ccLagPts,osi*.001,'intv',0);
    hold on
    for i=bpix
      if length(diag(r(i).gaeCCMn))>1
        sph=subplot('position',[xpos(lct) ypos(row) xlen ylen]);
        lct=lct+1;
        srfc=repmat(nan,length(cci),nAllLFPCh);
        % the lines below take into account non-analyzed channels (which have
        % single nans in the corresponding cells)
        for g=1:nAllLFPCh
          ix1=min(g,LFPpcInd2);
          ix2=max(g,LFPpcInd2);
          if length(r(i).gaeCCMn{ix1,ix2})>1
            srfc(1:length(cci),g)=r(i).gaeCCMn{ix1,ix2}(cci);
            % symmetry
            if g>LFPpcInd2, srfc(:,g)=flipud(srfc(:,g)); end
          else
            srfc(1:length(cci),g)=nan;
          end
        end
        % --- contour plot
        [c,cph]=contourf(1:nAllLFPCh,ccp,srfc,nLev);
        clear c;
        axis tight
        set(cph(:),'linestyle','none');
        cph=gca;
        % cbh=colorbar(cbPos);
        % use the contour plot's x axis as electrode # indicator
        set(cph,'xaxisloc','top','yaxisloc','right','ytick',[],...
          'xlim',[1 nAllLFPCh],'xtick',[1:nAllLFPCh],'xticklabel',WP.xtl2);
        % x axis limits must be fixed for aligned overlay with surface plot
        ax2=axes('position',get(cph,'position'),'color','none');
        hold on;
        set(gca,'ylim',[ccp(1) ccp(end)],'xlim',[WP.elx(1) WP.elx(end)],...
          'xtick',WP.elx,'xticklabel',WP.xtl1);
        axis manual;
        pcol=AP.segmentType{i,3}; psym=AP.segmentType{i,4};
        tmpMn=cat(1,r(i).gaeCCPeakTMn{1:LFPpcInd2,LFPpcInd2});
        % symmetry - don't forget to invert
        tmpMn=[tmpMn; -1*cat(1, r(i).gaeCCPeakTMn{LFPpcInd2,LFPpcInd2+1:end})];
        tmpStd=cat(1, r(i).gaeCCPeakTStd{1:LFPpcInd2,LFPpcInd2});
        tmpStd=[tmpStd; cat(1, r(i).gaeCCPeakTStd{LFPpcInd2,LFPpcInd2+1:end})];
        ph=errorbar(WP.elx,tmpMn,tmpStd,[psym '-']);
        set(ph,'color',pcol);
        if lct==2, ylabel('lag (ms)'); end
        grid on;
        % line to indicate principal electrode
        lh=line(WP.elx(LFPpcInd2)*[1 1]',get(gca,'ylim'));
        set(lh,'color',pcCol,'linestyle',':');
        % push to background
        set(gca,'children',circshift(get(gca,'children'),-1))
        title(['XC{\gamma}Env,{\gamma}Env(ref);' r(i).segmentType]);
      end
    end
  end

  if isfield(r,'thLoeCCMn')
    % II. thetaLo env (range to plot: the lesser of [+/- 1 slow theta period, AP.ccLagPts])
    row=2; lct=1;
    ccw=[-1 1]*min(AP.ccLagPts,cont2discrete(1*1000/AP.theta(1),osi*.001,'intv',1));
    cci=AP.ccLagPts+[ccw(1):ccw(2)]+1;
    ccp=discrete2cont(cci-AP.ccLagPts,osi*.001,'intv',0);
    hold on
    for i=bpix
      if length(diag(r(i).thLoeCCMn))>1
        sph=subplot('position',[xpos(lct) ypos(row) xlen ylen]);
        lct=lct+1;
        srfc=repmat(nan,length(cci),nAllLFPCh);
        % the lines below take into account non-analyzed channels (which have
        % single nans in the corresponding cells)
        for g=1:nAllLFPCh
          ix1=min(g,LFPpcInd2);
          ix2=max(g,LFPpcInd2);
          if length(r(i).thLoeCCMn{ix1,ix2})>1
            srfc(1:length(cci),g)=r(i).thLoeCCMn{ix1,ix2}(cci);
            % symmetry
            if g>LFPpcInd2, srfc(:,g)=flipud(srfc(:,g)); end
          else
            srfc(1:length(cci),g)=nan;
          end
        end
        % --- contour plot
        [c,cph]=contourf(1:nAllLFPCh,ccp,srfc,nLev);
        clear c;
        axis tight
        set(cph(:),'linestyle','none');
        cph=gca;
        % cbh=colorbar(cbPos);
        % use the contour plot's x axis as electrode # indicator
        set(cph,'xaxisloc','top','yaxisloc','right','ytick',[],...
          'xlim',[1 nAllLFPCh],'xtick',[1:nAllLFPCh],'xticklabel',WP.xtl2);
        % x axis limits must be fixed for aligned overlay with surface plot
        ax2=axes('position',get(cph,'position'),'color','none');
        hold on;
        set(gca,'ylim',[ccp(1) ccp(end)],'xlim',[WP.elx(1) WP.elx(end)],...
          'xtick',WP.elx,'xticklabel',WP.xtl1);
        axis manual;
        pcol=AP.segmentType{i,3}; psym=AP.segmentType{i,4};
        tmpMn=cat(1,r(i).thLoeCCPeakTMn{1:LFPpcInd2,LFPpcInd2});
        % symmetry - don't forget to invert
        tmpMn=[tmpMn; -1*cat(1, r(i).thLoeCCPeakTMn{LFPpcInd2,LFPpcInd2+1:end})];
        tmpStd=cat(1, r(i).thLoeCCPeakTStd{1:LFPpcInd2,LFPpcInd2});
        tmpStd=[tmpStd; cat(1, r(i).thLoeCCPeakTStd{LFPpcInd2,LFPpcInd2+1:end})];
        ph=errorbar(WP.elx,tmpMn,tmpStd,[psym '-']);
        set(ph,'color',pcol);
        if lct==2, ylabel('lag (ms)'); end
        grid on;
        % line to indicate principal electrode
        lh=line(WP.elx(LFPpcInd2)*[1 1]',get(gca,'ylim'));
        set(lh,'color',pcCol,'linestyle',':');
        % push to background
        set(gca,'children',circshift(get(gca,'children'),-1))
        title(['XC{\theta}LoEnv,{\theta}LoEnv(ref);' r(i).segmentType]);
      end
    end
  end


  if isfield(r,'thHieCCMn')
    % III. thetaHi env (range to plot: the lesser of [+/- 1 slow theta period, AP.ccLagPts])
    row=3; lct=1;
    ccw=[-1 1]*min(AP.ccLagPts,cont2discrete(1*1000/AP.theta(1),osi*.001,'intv',1));
    cci=AP.ccLagPts+[ccw(1):ccw(2)]+1;
    ccp=discrete2cont(cci-AP.ccLagPts,osi*.001,'intv',0);
    hold on
    for i=bpix
      if length(diag(r(i).thHieCCMn))>1
        sph=subplot('position',[xpos(lct) ypos(row) xlen ylen]);
        lct=lct+1;
        srfc=repmat(nan,length(cci),nAllLFPCh);
        % the lines below take into account non-analyzed channels (which have
        % single nans in the corresponding cells)
        for g=1:nAllLFPCh
          ix1=min(g,LFPpcInd2);
          ix2=max(g,LFPpcInd2);
          if length(r(i).thHieCCMn{ix1,ix2})>1
            srfc(1:length(cci),g)=r(i).thHieCCMn{ix1,ix2}(cci);
            % symmetry
            if g>LFPpcInd2, srfc(:,g)=flipud(srfc(:,g)); end
          else
            srfc(1:length(cci),g)=nan;
          end
        end
        % --- contour plot
        [c,cph]=contourf(1:nAllLFPCh,ccp,srfc,nLev);
        clear c;
        axis tight
        set(cph(:),'linestyle','none');
        cph=gca;
        % cbh=colorbar(cbPos);
        % use the contour plot's x axis as electrode # indicator
        set(cph,'xaxisloc','top','yaxisloc','right','ytick',[],...
          'xlim',[1 nAllLFPCh],'xtick',[1:nAllLFPCh],'xticklabel',WP.xtl2);
        % x axis limits must be fixed for aligned overlay with surface plot
        ax2=axes('position',get(cph,'position'),'color','none');
        hold on;
        set(gca,'ylim',[ccp(1) ccp(end)],'xlim',[WP.elx(1) WP.elx(end)],...
          'xtick',WP.elx,'xticklabel',WP.xtl1);
        axis manual;
        pcol=AP.segmentType{i,3}; psym=AP.segmentType{i,4};
        tmpMn=cat(1,r(i).thHieCCPeakTMn{1:LFPpcInd2,LFPpcInd2});
        % symmetry - don't forget to invert
        tmpMn=[tmpMn; -1*cat(1, r(i).thHieCCPeakTMn{LFPpcInd2,LFPpcInd2+1:end})];
        tmpStd=cat(1, r(i).thHieCCPeakTStd{1:LFPpcInd2,LFPpcInd2});
        tmpStd=[tmpStd; cat(1, r(i).thHieCCPeakTStd{LFPpcInd2,LFPpcInd2+1:end})];
        ph=errorbar(WP.elx,tmpMn,tmpStd,[psym '-']);
        set(ph,'color',pcol);
        if lct==2, ylabel('lag (ms)'); end
        grid on;
        % line to indicate principal electrode
        lh=line(WP.elx(LFPpcInd2)*[1 1]',get(gca,'ylim'));
        set(lh,'color',pcCol,'linestyle',':');
        % push to background
        set(gca,'children',circshift(get(gca,'children'),-1))
        title(['XC{\theta}HiEnv,{\theta}HiEnv(ref);' r(i).segmentType]);
      end
    end
  end

  if length(get(gcf,'children'))<=1
    % this prevents empty figures from remaining on screen and getting printed
    delete(gcf);
  else
    % print this summary figure?
    if ~isempty(AP.printas{1}),
      for i=1:length(AP.printas)
        pa=AP.printas{i};
        if strfind(pa,'ps'), ext='.ps';
        elseif strfind(pa,'jpeg'), ext='.jpg';
        else ext='';
        end
        print(pa,[figName '_ccContour2' ext]);
      end
    end
  end
end



% -------------------------------------------------------------------------------
%               ***  FIGURE V: images of correlation matrices ***
% -------------------------------------------------------------------------------
fh5=mkfig('CCmatImg');
orient tall
labelscale('fontSz',7,'scaleFac',1.0,'lineW',1.0,'markSz',4);
% set up compound colormap
nColors=128;
colormap([jet(nColors); coma('redblue','ncols',nColors)]);
strmTypes={'theta','gamma','thetaLoEnv','thetaHiEnv','gammaEnv','delta'};
% # of rows = # of CC types
nRows=length(strmTypes);
% # of cols = # of behavior types with results * 2
nCols=2*length(bpix);
xmarg=.04;
ymarg=.035;
xpos=(0:nCols-1)*((1-2*xmarg)/nCols)+2*xmarg;
xlen=(1-2*xmarg)/nCols-2*xmarg;
% leave additional space (.05)for title subplot
ypos=fliplr((0:nRows-1)*((1-.05-2*ymarg)/nRows)+2*ymarg);
ylen=(1-.05-2*ymarg)/nRows-2*ymarg;
% title
subplot('position',[xmarg .95 1-2*xmarg .02]);
set(gca,'xlim',[0 1],'ylim',[0 1]);
axis off
th=text(.5,.5,figTitle,'fontsize',12,'fontweight','bold','horizontalalignment','center');
ccprop={'ampl','lag'};
cblabFstr={'%1.1f','%3.0f'};
ccMatTemplate=repmat(nan,[nAllLFPCh nAllLFPCh]);

for sti=1:nRows
  strmType=strmTypes{sti};
  switch strmType
    case 'delta'
      maxLag=250;
      minAmp=.0;
      % this stream's short form used in field of struct r
      STshort='de';
    case 'theta'
      maxLag=100;
      minAmp=.0;
      STshort='th';
    case 'thetaLoEnv'
      maxLag=180;
      minAmp=.0;
      % this stream's short form used in field of struct r
      STshort='thLoe';
    case 'thetaHiEnv'
      maxLag=180;
      minAmp=.0;
      % this stream's short form used in field of struct r
      STshort='thHie';
    case 'gamma'
      maxLag=15;
      minAmp=.0;
      % this stream's short form used in field of struct r
      STshort='ga';
    case 'gammaEnv'
      maxLag=15;
      minAmp=.0;
      % this stream's short form used in field of struct r
      STshort='gae';
    otherwise
      error('illegal streamType in CC matrix plot');
  end

  % now the big scaling: since a colormap applies to a whole figure (and not a
  % single axis) we have to use a compound color map and scale the data such that
  % CC lags and amplitudes are converted to sets of nonoverlapping indices
  % 1. amplitudes: map [0 1] to [1 nColors]
  cFac(1)=nColors-1;
  cOffs(1)=1;
  % 2. lags: map [-maxLag maxLag] to [nColors+1 2*nColors]
  cFac(2)=(nColors-1)/(2*maxLag);
  cOffs(2)=cFac(1)+cOffs(1)+nColors/2;
  ict=0;
  row=sti;
  if isfield(r,[STshort 'CCMn'])
    for i=bpix
      ict=ict+1;
      eval(['cm1=diag(r(i).' STshort 'CCMn);']);
      if length(cm1)>1
        eval(['cm1=r(i).' STshort 'CCPeakMn;']);
        eval(['cm2=r(i).' STshort 'CCPeakTMn;']);
        ccLag=ccMatTemplate;
        ccAmp=ccMatTemplate;
        % transfer values from cell array in array - has to be done elementwise because
        % there are empty cells
        for g=1:nAllLFPCh^2
          if ~isempty(cm1{g})
            ccAmp(g)=cm1{g};
            ccLag(g)=cm2{g};
          end
        end
        % 1. amplitudes
        % CC peak values should be within [0 1], but to prevent indexing errors make
        % sure that is really so:
        ccAmp(ccAmp>1)=1;
        ccAmp(ccAmp<0)=0;
        ccAmp=round(cOffs(1)+cFac(1)*ccAmp);
        % 2. lags:
        % same story: deal with outliers (set to nan here)
        ccLag(abs(ccLag)>maxLag)=nan;
        ccLag=round(cOffs(2)+cFac(2)*ccLag);
        % ---
        for ii=1:2
          sph=subplot('position',[xpos(ict+(ii-1)*2) ypos(row) xlen ylen]);
          if ii==1
            ih=image(ccAmp);
          else
            ih=image(ccLag);
          end
          % set(ih,'CDataMapping','scaled')
          set(ih,'CDataMapping','direct')
          cph=gca;
          set(cph,'xaxisloc','top',...
            'ylim',[.5 nAllLFPCh+.5],'ytick',[1:nAllLFPCh],'yticklabel',WP.xtl2,...
            'xlim',[.5 nAllLFPCh+.5],'xtick',[1:nAllLFPCh],'xticklabel',WP.xtl2);
          axis square
          % lines to indicate principal electrode
          lh=line(LFPpcInd2*[1 1]',[.5 LFPpcInd2]);
          set(lh,'color',[.6 .6 .6],'linestyle','-');
          lh=line([LFPpcInd2 nAllLFPCh+.5],LFPpcInd2*[1 1]');
          set(lh,'color',[.6 .6 .6],'linestyle','-');
          cbh=colorbar(cbPos);
          set(cbh,'ylim',[(ii-1)*nColors  ii*nColors]);
          % set ticks/labels at/to meaningful positions/values
          set(cbh,'ytick',[ii-1 ii-.5  ii]*nColors,...
            'yticklabel',num2str(((([ii-1 ii-.5  ii]*nColors)-cOffs(ii))/cFac(ii))',cblabFstr{ii}));
          % inflate plot
          rexy('ax',cph,'xfac',1.3,'yfac',1.3);
          title([r(i).segmentType ';' strmType ';' ccprop{ii}]);
        end
      end
    end
  end
end

if length(get(gcf,'children'))<=1
  % this prevents empty figures from remaining on screen and getting printed
  delete(gcf);
else
  % print this summary figure?
  if ~isempty(AP.printas{1}),
    for i=1:length(AP.printas)
      pa=AP.printas{i};
      if strfind(pa,'ps'), ext='.ps';
      elseif strfind(pa,'jpeg'), ext='.jpg';
      else ext='';
      end
      print(pa,[figName '_ccMatImg' ext]);
    end
  end
end


% -------------------------------------------------------------------------------
%    ***  FIGURE VI: comodulograms ***
% -------------------------------------------------------------------------------

fh6=mkfig('fComod');
orient tall
labelscale('fontSz',6,'scaleFac',1.0,'lineW',1.0,'markSz',4);
colormap(flipud(coma('amberturquois','ncols',255)));
% assume that most correlation coeffs are in that range
ccLim=[-.5 .5];
% # of rows = max # of channels /2
nRows=ceil(length(AP.allLFPIdx)/2);
% # of cols = # of behavior types with results *2
nCols=length(bpix)*2;
xmarg=.04;
ymarg=.035;
xpos=(0:nCols-1)*((1-2*xmarg)/nCols)+2*xmarg;
xlen=(1-2*xmarg)/nCols-2*xmarg;
% leave additional space (.05)for title subplot
ypos=fliplr((0:nRows-1)*((1-.05-2*ymarg)/nRows)+2*ymarg);
ylen=(1-.05-2*ymarg)/nRows-2*ymarg;
% title
subplot('position',[xmarg .95 1-2*xmarg .02]);
set(gca,'xlim',[0 1],'ylim',[0 1]);
axis off
th=text(.5,.5,figTitle,'fontsize',12,'fontweight','bold','horizontalalignment','center');

if isfield(r,'fComod')
  lct=0;
  for i=bpix
    if ~isempty(r(i).fComod)
      lct=lct+1;
      for tmpi=1:nLFPCh
        row=mod(AP.LFPccInd(tmpi)-1,nRows)+1;
        if AP.LFPccInd(tmpi)==nRows+1, lct=lct+1; end
        sph=subplot('position',[xpos(lct) ypos(row) xlen ylen]);
        cph=imagesc(r(1).comF,r(1).comF,r(i).fComod(:,:,tmpi),ccLim);
        axis square
        % inflate plot
        rexy('ax',gca,'xfac',1.1,'yfac',1.1);
        set(gca,'xtick',[10:20:r(1).comF(end)],'ytick',[10:20:r(1).comF(end)]);
        grid on
        title(AP.rawChAnNm{tmpi});
      end
    end
  end
end

% print this summary figure?
if ~isempty(AP.printas{1}),
  for i=1:length(AP.printas)
    pa=AP.printas{i};
    if strfind(pa,'ps'), ext='.ps';
    elseif strfind(pa,'jpeg'), ext='.jpg';
    else ext='';
    end
    print(pa,[figName '_fComod' ext]);
  end
end

% -------------------------------------------------------------------------------
%    ***  FIGURE VIb: prob values of comodulograms ***
% -------------------------------------------------------------------------------

fh6=mkfig('fComodP');
orient tall
labelscale('fontSz',6,'scaleFac',1.0,'lineW',1.0,'markSz',4);
colormap(gray);
pLim=[0 .05];
% # of rows = max # of channels /2
nRows=ceil(length(AP.allLFPIdx)/2);
% # of cols = # of behavior types with results *2
nCols=length(bpix)*2;
xmarg=.04;
ymarg=.035;
xpos=(0:nCols-1)*((1-2*xmarg)/nCols)+2*xmarg;
xlen=(1-2*xmarg)/nCols-2*xmarg;
% leave additional space (.05)for title subplot
ypos=fliplr((0:nRows-1)*((1-.05-2*ymarg)/nRows)+2*ymarg);
ylen=(1-.05-2*ymarg)/nRows-2*ymarg;
% title
subplot('position',[xmarg .95 1-2*xmarg .02]);
set(gca,'xlim',[0 1],'ylim',[0 1]);
axis off
th=text(.5,.5,figTitle,'fontsize',12,'fontweight','bold','horizontalalignment','center');

if isfield(r,'fComodP')
  lct=0;
  for i=bpix
    if ~isempty(r(i).fComodP)
      lct=lct+1;
      for tmpi=1:nLFPCh
        row=mod(AP.LFPccInd(tmpi)-1,nRows)+1;
        if AP.LFPccInd(tmpi)==nRows+1, lct=lct+1; end
        sph=subplot('position',[xpos(lct) ypos(row) xlen ylen]);
        cph=imagesc(r(1).comF,r(1).comF,r(i).fComodP(:,:,tmpi),pLim);
        axis square
        % inflate plot
        rexy('ax',gca,'xfac',1.1,'yfac',1.1);
        set(gca,'xtick',[10:20:r(1).comF(end)],'ytick',[10:20:r(1).comF(end)]);
        grid on
        title(AP.rawChAnNm{tmpi});
      end
    end
  end
end

% print this summary figure?
if ~isempty(AP.printas{1}),
  for i=1:length(AP.printas)
    pa=AP.printas{i};
    if strfind(pa,'ps'), ext='.ps';
    elseif strfind(pa,'jpeg'), ext='.jpg';
    else ext='';
    end
    print(pa,[figName '_fComodP' ext]);
  end
end

% % -------------------------------------------------------------------------------
% %    ***  FIGURE VIc: freq band-specific means of comodulograms ***
% % -------------------------------------------------------------------------------
% 
% fh6c=mkfig('fComodSum');
% orient tall
% labelscale('fontSz',6,'scaleFac',1.0,'lineW',1.0,'markSz',4);
% colormap(flipud(coma('amberturquois','ncols',255)));
% % assume that most correlation coeffs are in that range
% ccLim=[-.5 .5];
% % # of rows = max # of channels /2
% nRows=ceil(length(AP.allLFPIdx)/2);
% % # of cols = # of behavior types with results *2
% nCols=length(bpix)*2;
% xmarg=.04;
% ymarg=.035;
% xpos=(0:nCols-1)*((1-2*xmarg)/nCols)+2*xmarg;
% xlen=(1-2*xmarg)/nCols-2*xmarg;
% % leave additional space (.05)for title subplot
% ypos=fliplr((0:nRows-1)*((1-.05-2*ymarg)/nRows)+2*ymarg);
% ylen=(1-.05-2*ymarg)/nRows-2*ymarg;
% 
% 
% if isfield(r,'fComod')
%   fBand=[AP.delta; AP.theta; AP.beta; AP.gamma];
%   fbLabel={'\delta','\theta','\beta','\gamma'};
%   for fbi=1:size(fBand,1)
%     fBandFIx{fbi}=find(r(1).comF>=fBand(fbi,1) & r(1).comF<=fBand(fbi,2));
%   end
%   fbMat=repmat(0,length(fBandFIx)*[1 1]);
% 
%   lct=0;
%   for i=bpix
%     if ~isempty(r(i).fComod)
%       lct=lct+1;
%       for tmpi=1:nLFPCh
%         row=mod(AP.LFPccInd(tmpi)-1,nRows)+1;
%         if AP.LFPccInd(tmpi)==nRows+1, lct=lct+1; end
%         sph=subplot('position',[xpos(lct) ypos(row) xlen ylen]);
%         hold on
%         
%         for fbi1=1:length(fBandFIx)
%           for fbi2=fbi1:length(fBandFIx)
%             tmpMat=r(i).fComod(fBandFIx{fbi1},fBandFIx{fbi2},tmpi);
%             fbMat(fbi1,fbi2)=mean(tmpMat(:));
%           end
%           th=text(fbi1-.5,.5,fbLabel{fbi1});
%           th=text(.5,fbi1-.5,fbLabel{fbi1});
%         end
%         
%         imagesc(fbMat,ccLim);
%         set(gca,'ydir','reverse');
%         axis square
%         % inflate plot
%         rexy('ax',gca,'xfac',1.5,'yfac',1.5);
%         set(gca,'xtick',[],'ytick',[]);
% %         tmpMn=mean(r(i).fComod(:,:,row),1);
% %         ph=plot(r(1).comF,tmpMn,'k-');
% %         ph=plot(r(1).comF,tmpMn+std(r(i).fComod(:,:,row),0,1),'k:');        
% %         % same, but only significant pixels
% %         tmpFix=(r(i).fComodP(:,:,row)<=.05);
% %         for g=1:length(r(1).comF)
% %           tmpMn(g)=mean(r(i).fComod(find(tmpFix(:,g)),g,row));
% %           tmpStd(g)=std(r(i).fComod(find(tmpFix(:,g)),g,row));          
% %         end
% %         ph=plot(r(1).comF,tmpMn,'m-');
% %         ph=plot(r(1).comF,tmpMn+std(r(i).fComod(:,:,row),0,1),'m:');        
% %         niceyax;
% %         line(get(gca,'xlim'),[0 0]);
% % 
% %         % set(gca,'ylim',ccLim);
% %         % inflate plot
% %         rexy('ax',gca,'xfac',1.2,'yfac',1.5);
% %         set(gca,'xtick',[10:20:r(1).comF(end)],'ytick',[ccLim(1):.1:ccLim(2)]);
%         title(AP.rawChAnNm{tmpi});
%       end
%     end
%   end
% end
% % subpax(fh6b);
% % % title - has to come last here because of subpax
% subplot('position',[xmarg .95 1-2*xmarg .02]);
% set(gca,'xlim',[0 1],'ylim',[0 1]);
% axis off
% th=text(.5,.5,figTitle,'fontsize',12,'fontweight','bold','horizontalalignment','center');
% 
% % print this summary figure?
% if ~isempty(AP.printas{1}),
%   for i=1:length(AP.printas)
%     pa=AP.printas{i};
%     if strfind(pa,'ps'), ext='.ps';
%     elseif strfind(pa,'jpeg'), ext='.jpg';
%     else ext='';
%     end
%     print(pa,[figName '_fComodSum' ext]);
%   end
% end

% list not complete
clear tmp* ph sph cc_contourP spec_contourP temnpo pcol frix fifi titulo th tp lh xl yl pa ext cph cbh row lct g ix* nLev nRows nCols narrThF ppos fBand

warning on

% ---------------------------------- crap ---------------------------
