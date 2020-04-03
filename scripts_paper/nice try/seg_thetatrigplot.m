% This routine generates plots of theta, gamma and raw data triggered to
% individual peaks of theta. It is a minor modification of
% seg_thetatrigstreams.m; therefore, it must be invoked in the exact same
% way
AP.printas={'-dpsc2'};
figdir='d:\projects\rmouse\paper_atropine\rawFig\';
% figdir='';
figName=[mfilename '_' DS.abfFn(end)]; 

thetaPowCo=.10;
thetaPeakCo=.10;
thetaDeadT=50;
% peri-peak interval to cut out (ms)
excIntv=[-175 175];
excIntvPts=cont2discrete(excIntv,osi*.001,'intv',1);
eIx=excIntvPts(1):excIntvPts(2);
% borders to cut off from excerpts (filter) in ms
border=50;
killPts=cont2discrete(border,osi*.001);
killPts=[1:killPts length(eIx)-killPts+1:length(eIx)];
% max # of theta peaks to analyze
maxNThetaPeak=5000;
% max # of traces to plot
maxNThetaPeakPlot=200;


% choose channels to analyze: slm (principal) and oriens (i.e. 500 um dorsad)
dDist=.5;
[doff,ccRefChInd]=min(abs(WP.elx(AP.LFPccInd)+dDist));
ccRefChIdx=AP.LFPIdx(ccRefChInd);
elPos=WP.elx([ccRefChIdx AP.pcIdx]);
dDist=diff(elPos);
nChannel=2;

locChIdx=[ccRefChIdx AP.LFPpcInd2];
locChInd=[ccRefChInd AP.LFPpcInd1];

% graphics
fh=mkfig('thetaTrigRaw');
orient portrait
labelscale('fontSz',10,'scaleFac',1,'lineW',.75,'markSz',4);

% ** exploring only **
for i=explV
  if ~isempty(r(i).iPts)
    disp([r(i).segmentType '..']);

    % load theta peaks (time is ms)
    load(WP.thetaPFn,'negPA','negPT','posPA','posPT');
    % subplot column index
    spci=0;
    for chaI=1:2
      chIdx=locChIdx(chaI);
      chInd=locChInd(chaI);      
      
      disp([r(i).segmentType ': theta peak-triggered raw plots, ch ' rawCh(chInd).nm ' ..']);
      spci=spci+1;
      % identify segments with theta power above threshold
      [tmpCo]=cumh((r(i).rawThPE{chIdx,chIdx})',.01,'p',thetaPowCo);
      segIx=find(r(i).rawThPE{chIdx,chIdx}>=tmpCo);
      nSeg=length(segIx);
      
      if chaI==1
        % dorsal sites: trigger to positive peaks
        peakA=posPA{chIdx};
        peakT=posPT{chIdx};
      else
        % close to fissure: trigger to negative peaks
        peakA=-1*negPA{chIdx};
        peakT=negPT{chIdx};
      end
        
      % identify theta peaks above threshold
      [tmpCo]=cumh(peakA,.01,'p',thetaPeakCo);

      % time stamp list of potentially interesting peaks in pts
      thPeakIx=(peakA>=tmpCo);
      thPeakT=peakT(thPeakIx);
      % purge double peaks
      thPeakT=tsldeadt(thPeakT,thetaDeadT);
      % convert to pts
      thPeakT=cont2discrete(thPeakT,osi*.001);
      % remove peaks too close to borders
      thPeakT(thPeakT+eIx(1)<1)=[];
      thPeakT(thPeakT+eIx(end)>r(i).iPts(segIx(end),2))=[];
      
      % all remaining peaks in segments
      thPeakIx=zeros(size(thPeakT));
      for g=1:nSeg
        thPeakIx=thPeakIx | (thPeakT>=r(i).iPts(segIx(g),1) & thPeakT<=r(i).iPts(segIx(g),2));
      end
      thPeakT=thPeakT(thPeakIx);
      % reduce to max number of peaks to be analyzed
      step=max(1,length(thPeakT)/maxNThetaPeak);
      thPeakT=thPeakT(round(1:step:length(thPeakT)));
      nThPeak=length(thPeakT);

      % indices to peaks to be plotted
      step=max(1,nThPeak/maxNThetaPeakPlot);
      plotIx=round(1:step:nThPeak);
      
      % load raw data: beginning to last point of last identified segment
      if exist([dpath '\' abfFn '.mat'],'file')
        rawD=matDload([dpath '\' abfFn '.mat'],'start',0,'stop',discrete2cont(r(i).iPts(segIx(end),2),osi*1e-6,'intv',1),...
          'channels',{rawCh(chInd).nm},'verbose',verbose);
      else
        rawD=abfload([dpath '\' abfFn '.abf'],'start',0,'stop',discrete2cont(r(i).iPts(segIx(end),2),osi*1e-6,'intv',1),...
          'channels',{rawCh(chInd).nm},'verbose',verbose);
      end
      if DS.rawSignalInverted
        rawD=-1*rawD;
      end
      % Due to the fact that abfload computes discrete time from continuous time a little
      % differently than is done in rmouse_xx, in some cases a single data point is
      % amiss. This is bad, and calls for a revision of abfload. In the
      % meantime, here's the workaround:
      tmppd=r(i).iPts(segIx(end),2)-length(rawD);
      if tmppd,
        if tmppd==1
          disp('missed a point');
          rawD(end+1)=rawD(end);
        elseif tmppd== -1
          disp('one point too much');
        else
          error('abfload or matDload handles time information poorly');
        end
      end
      % preallocate
      excD=repmat(nan,length(eIx),nThPeak-2);
      for g=1:nThPeak
        excD(:,g)=rawD(eIx+thPeakT(g));
      end
      % filter:
      % 'raw' = 4-90
      excD=bafi(excD,si,[AP.thetaCFreq(1) AP.gammaCFreq(2)],'rs',AP.rs);
      % gamma as usual 
      gaD=bafi(excD,si,AP.gammaCFreq,'rs',AP.rs);
      % 'fast': everything above upper theta and below upper gamma
      fastD=bafi(excD,si,[AP.thetaCFreq(2) AP.gammaCFreq(2)],'rs',AP.rs);
      % cut
      excD(killPts,:)=[];
      gaD(killPts,:)=[];
      fastD(killPts,:)=[];
      % offsets
%       tmpd=excD(:,plotIx);
%       co1=cumh(tmpd(:),.001,'p',[.001]);
%       tmpd=fastD(:,plotIx);
%       co2=cumh(tmpd(:),.001,'p',[.001 .999]);
%       tmpd=gaD(:,plotIx);
%       co3=cumh(tmpd(:),.001,'p',[.999]);
%       offs=cumsum([0 co1-co2(2)  co2(1)-co3])

      if chaI==1
        offs=[0  -1.4   -2.5];
        ylim=[-3.0 .9];
      else
        offs=[0  -3.8   -6.4];
        ylim=[-7.2 1.6];
      end


%       tmpd=beD(:,plotIx);
%       co2=cumh(tmpd(:),.001,'p',[.001 .999]);
%       tmpd=gaD(:,plotIx);
%       co3=cumh(tmpd(:),.001,'p',[.999]);
%       offs=cumsum([0 co1-co2(2)  co2(1)-co3]);
      
      % stack plots in one row
      subplot(nChannel,1,spci), hold on
      rexy('ax','gca','xfac',.22,'yfac',1);
      
      ph=plot(excD(:,plotIx)+offs(1));
      set(ph,'color',[.6 .6 .6]);
      
      ph=plot(fastD(:,plotIx)+offs(2));
      set(ph,'color',[.6 .6 .6]);
      
      ph=plot(gaD(:,plotIx)+offs(3));
      set(ph,'color',[.6 .6 .6]);

      ph=plot(mean(excD,2)+offs(1),'k');      
      set(ph,'linewidth',1.6);
      ph=plot(mean(fastD,2)+offs(2),'k');      
      set(ph,'linewidth',1.6);
      ph=plot(mean(gaD,2)+offs(3),'k');      
      set(ph,'linewidth',1.6);
      
      ph=plot(mean(excD-fastD,2)+offs(1),'c');      
      set(ph,'linewidth',1);

      
      title(['rec pos ' num2str(elPos(spci)*1000,'%4.0f') ' um']);

      % axis tight;
      set(gca,'ylim',ylim);
      axis off;
      
      utscaleb3(osi);
      
    end % for:channels
  end
end

drawnow
if ~isempty(AP.printas{1}),
  for i=1:length(AP.printas)
    pa=AP.printas{i};
    if strfind(pa,'ps'), ext='.ps';
    elseif strfind(pa,'jpeg'), ext='.jpg';
    else ext='';
    end
    print(pa,[figdir figName ext]);
  end
end
% reset graphics to default
labelscale('fontSz',7,'scaleFac',1.0,'lineW',1.0,'markSz',4);

