function rmouse_manualccpeak
% **** NOTE: this function will only run if AP.job contains 'seg_cc_ms2rad'
% (or any other 'seg_..' job) and AP.saveMode='append'!!

global DS AP WP r

% -------------- preliminaries ---------------------------

% --- prompt user for type of stream and behavior
disp('---------- manual detection of CC peaks..');
prompt={'Enter legal stream name (theta, gammaNarrow, etc.',...
  'Which behavior (immobile,exploring)?'};
name='Manual detection of CC peaks';
numlines=1;
default={'theta','exploring'};
answer=inputdlg(prompt,name,numlines,default);
strmType=answer{1};
behavInd=strmatch(answer{2},AP.segmentType(:,1),'exact');
% default entries into defExpLag
defaultLag=0;

% ---- timing, conversion discrete <-> continuous
switch strmType
  case 'delta'
    STshort='de';
    % detect peaks within this interval centered around t=0
    ccw=cont2discrete(AP.deccw,WP.osi*.001,'intv',1);
  case 'theta'
    STshort='th';
    ccw=cont2discrete(AP.thccw,WP.osi*.001,'intv',1);
  case 'thetaHiEnv'
    STshort='thHie';
    ccw=cont2discrete(AP.deccw,WP.osi*.001,'intv',1);
  case 'thetaLoEnv'
    STshort='thLoe';
    ccw=cont2discrete(AP.ccLagPts*[-.95 .95],WP.osi*.001,'intv',1);
  case {'gamma'}
    STshort='ga';
    ccw=cont2discrete(AP.gaccw,WP.osi*.001,'intv',1);
  case {'gammaEnv'}
    STshort='gae';
    ccw=cont2discrete(AP.gaccw,WP.osi*.001,'intv',1);
  case {'gammaNarrow'}
    STshort='gaNa';
    ccw=cont2discrete(AP.gaNaccw,WP.osi*.001,'intv',1);
  case {'gammaNarrowEnv'}
    STshort='gaNae';
    ccw=cont2discrete(AP.gaNaccw,WP.osi*.001,'intv',1);
  otherwise
    error('illegal streamType in manual CC peak determination');
end
% lag axis (in pts)
tmpAx=ccw(1):ccw(2);
% lag axis in ms
ax=discrete2cont(tmpAx+1,WP.osi*.001);
% index into excerpt of cc as defined in AP
cci=AP.ccLagPts+1+tmpAx;

% --- preallocation and other prep of variables
% array to hold determined lags
defExpLag=repmat(defaultLag,AP.nAllLFPCh,AP.nAllLFPCh)+tril(zeros(AP.nAllLFPCh)*nan,-1);

% generate index into cc matrix for upper triangular part going down columns 
for g=1:AP.nAllLFPCh
  chanPairInd{g}=[(1:g)' repmat(g,g,1) ];
end

% ------ mat file into which user-defined attractor lags (expLag) will be
% written
lagFn=[DS.abfFn '_userdef_' answer{2} '_' strmType 'cclag.mat'];

% ----- figure
labelscale('fontSz',10,'scaleFac',1.0,'lineW',1.0,'markSz',5);
tmpftag='CCLags';
fh_ccLag=findobj('tag',tmpftag);
if isempty(fh_ccLag)
  fh_ccLag=figure;
else
  figure(fh_ccLag);
end
tmpScrSz=get(0,'Screensize');
tmpScrSz(2)=tmpScrSz(4)*.2;
tmpScrSz(4)=tmpScrSz(4)*.7;
set(fh_ccLag,'position',tmpScrSz,'tag',tmpftag,'name','manual peak detection in CC' ,...
  'color',[0.8 0.8 0.8],'numbertitle','off');
% DONE button
h1 = uicontrol('Parent',fh_ccLag, ...
  'Units','normalized', ...
  'HorizontalAlignment','left', ...    
  'Position',[.8 .05 .1 .1], ...
  'FontSize',12, ...
  'Fontweight','bold',...
  'style','pushbutton',...
  'String','DONE', ...
  'ListboxTop',0, ...
  'Tag','doneBttn', ...
  'callback', @donefunc);  

% nested func: callback of DONE button (save & clean up)
  function donefunc(src,evdat)
    save(lagFn,'defExpLag');
    % reset graphics settings
    labelscale('fontSz',7,'scaleFac',1.0,'lineW',1.0,'markSz',4);
    if ishandle(fh_ccLag),
      delete(fh_ccLag);
    end
  end

% nested func: callback of peaks
  function togglePointSelect(src,evdat)
    uData=get(src,'userdata');
    if ~uData(4)
      % select
      uData(4)=1;
      set(src,'markerfacecolor',selectColor);
      defExpLag(uData(1),uData(2))=uData(3);
    else
      % unselect
      uData(4)=0;
      set(src,'markerfacecolor','none');
      defExpLag(uData(1),uData(2))=defaultLag;
    end
    set(src,'userdata',uData);
  end

% ------ further graphics settings
% y offset of cc functions
yOffs=1.5;
% markerfacecolor for selected peaks
selectColor=[1 0 0];

% ------ do it!
% (philosophy: plot in real time units (ms) but compute peaks on basis of
% ticks and then convert to ms)
cccMn=eval(['r(behavInd).' STshort 'CCMn;']);
% *** loop over diagonals ***
for ki=1:AP.nAllLFPCh
  ccix=chanPairInd{ki};
  nPairs=size(chanPairInd{ki},1);
  subplot(1,AP.nAllLFPCh,ki)
  hold on
  set(gca,'ytick',[],'color',[.8 .8 .8]);
  rexy('xfac',1.05,'yfac',.75);
  % leave space (in y dir) for peaks of 1 and assume an offset of 1.5
  % between traces
  axis([ax([1 end]) -AP.nAllLFPCh*1.5 1.5]);
  axis manual
  grid on
  title(DS.rawCh{AP.allLFPIdx(ki),1},'interpreter','none')
  % loop over elements in diagonals
  for j=1:nPairs
    ccix1=chanPairInd{ki}(j,1);
    ccix2=chanPairInd{ki}(j,2);
    if all(~isnan(cccMn{ccix1,ccix2}))
      % plot
      ph=plot(ax,cccMn{ccix1,ccix2}(cci)-(j-1)*yOffs,'k-');
      if ccix2==AP.LFPpcInd2
        % plot principal channel in thick trace 
        set(ph,'linewidth',2);
      end
      tmpr=evdeal(cccMn{ccix1,ccix2}(cci),'idx','allpeaks');
      % local var & subtract time offset in pts
      posPeakT=tmpr.posPeakT{1}+tmpAx(1);
      % *** lags in real time units ***
      posPeakTMs=discrete2cont(posPeakT,WP.osi*.001);
      for h=1:numel(posPeakT)
        ph=plot(posPeakTMs(h),tmpr.posPeak{1}(h)-(j-1)*yOffs,'bo');
        % *** userdata: row, column, lag in pts, selected or not
        set(ph,'userdata',[ccix1 ccix2 posPeakTMs(h) 0],'ButtondownFcn',@togglePointSelect)
      end
    end % for:j=1:nPairs
  end
end

waitfor(fh_ccLag);

end

