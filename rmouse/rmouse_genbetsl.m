function etsl=rmouse_genbetsl(bstart,bstop,b_delta,trigThresh)
% ** function etsl=rmouse_genbetsl(bstart,bstop,b_delta,trigThresh)
% generates extended time stamp list of behavior from raw data file
% containing trigger pulses

global DS AP WP 

% ---- preliminaries -----------------------
etslconst;

% the trigger levels in trigThresh MUST be correctly sorted at this point
negTrigInd=find(trigThresh(:,1)<0);
posTrigInd=find(trigThresh(:,1)>=0);
negTrigThresh=trigThresh(negTrigInd,1);
posTrigThresh=trigThresh(posTrigInd,1);
nTT=size(negTrigThresh,1);
pTT=size(posTrigThresh,1);

% ---- job proper -----------------------
[d,si]=abfload([DS.dpath '\' AP.bScoreFn '.abf'],'start',bstart*.001,'stop',bstop*.001,'channels',AP.rawChScoreNm);
% it is assumed that the trigger signals are NEVER inverted (unlike raw
% data, which may be inverted)
% downsample - sth on the order of 200 Hz ~ si=5000 us should be sufficient
sampF=round(5000/si);
d=d(1:sampF:end);
nsi=sampF*si;
% filter?
if ~isempty(AP.bsCFreq)
  d=hifi(d,nsi,AP.bsCFreq);
end
% this is, in ticks, the maximal duration we would expect a 'flank' of the trigger
% pulses to be
trigTauPts=round(AP.trigTau*1e3/nsi);
% negative triggers - etsl not yet sorted!
netsl=sspulse2tsl(d,negTrigThresh,nsi,-1,negTrigInd,trigTauPts);
% positive triggers - etsl not yet sorted!
petsl=sspulse2tsl(d,posTrigThresh,nsi,1,posTrigInd,trigTauPts);
if isempty(netsl) & isempty(petsl)
  msg={['data file: ' AP.resPath '\' DS.abfFn '; scoring file: ' AP.resPath '\' AP.bScoreFn ' - no behavioral trigger detected. Trigger thresholds may be off. Analysis stopped.']};
  warndlg(msg);
  disp(msg);
  return
end

% combine & sort
etsl=sortrows([netsl; petsl],1);

% ------ interlude: plot raw scoring data, excerpts of it &
%        trigger levels as defined in AP.behavType
% cut out triggers in interval slightly bigger than trigTau
ctout=tsl2exc(d,nsi,{etsl(:,etslc.tsCol)},'win',AP.trigTau*[-1.2  1.2]);
% --- figure
tmpftag='bscore';
fh_bscore=findobj('tag',tmpftag);
if isempty(fh_bscore), 
  fh_bscore=figure;
else
  figure(fh_bscore);
end
tmpScrSz=get(0,'Screensize');
tmpScrSz([2 4])=round([tmpScrSz(4)*.4  tmpScrSz(4)*.58]);
set(fh_bscore,'position',tmpScrSz,'tag',tmpftag,'name','behavioral scoring trigger pulses' ,...
  'color',[0.8 0.8 0.8],'numbertitle','off','menubar','none');
% - plot of raw scoring data & thresholds
subplot('position',[.03 .53 .77 .4]); hold on,
% plot lines first
tmpti=[negTrigInd; posTrigInd];
% time axis
tmpt=(1:length(d))*nsi*1e-6/60;
for i=1:length(tmpti)
  lh=line(tmpt([1 end]),trigThresh(i,1)*[1 1],'color',AP.behavType{tmpti(i),3},'linewidth',.5);
end
ph=plot(tmpt,d,'k');
set(ph,'linewidth',.5);
niceyax;
xlabel('time (min)');
ylabel('trigger level (mV)');
title('behavioral trigger pulses - overview (all from beginning of scoring to end of excerpt)','fontsize',12);
set(gca,'color',[0.85 0.85 0.85]);
box on;
yl=get(gca,'ylim');
% - all-points histogram & thresholds
[tmpn,tmpx]=hist(d,(0:.01:1)*diff(yl)+yl(1));
% fake entries so that zero entries appear on log scale
tmpn(tmpn==0)=.01;
subplot('position',[.82 .53 .15 .4]);
ph=contourbarh(tmpx,tmpn,'color','k','linewidth',1.0);
set(gca,'xscale','log','ylim',yl,'ytick',[],'yaxislocation','right','xlim',[.01 max(tmpn)*1.2]);
xl=get(gca,'xlim');
for i=1:length(tmpti)
  lh=line(xl,trigThresh(i,1)*[1 1],'color',AP.behavType{tmpti(i),3},'linewidth',.5);
end
xlabel('N');
ylabel('trigger level (mV)');
title('all-points histogram','fontsize',12);
set(gca,'color',[0.85 0.85 0.85]);
box on;
% - plot of trigger pulse excerpts
subplot('position',[.03 .08 .33 .37]); hold on,
ctout_t=cont2discrete([-1.2 1.2]*AP.trigTau,nsi*.001,'intv',1)+1;
plot((ctout_t(1):ctout_t(2))*nsi*.001,ctout,'o-');
niceyax;
xlabel('time (ms)');
ylabel('trigger level (mV)');
title('excerpts','fontsize',12);
set(gca,'color',[0.85 0.85 0.85]);

% UIs
hui0 = uicontrol('Parent',fh_bscore, ...
  'Units','normalized', ...
  'HorizontalAlignment','left', ...
  'Position',[.5 .35 .38 .07], ...
  'FontSize',12, ...
  'FontWeight','bold', ...
  'BackgroundColor',[0.8 0.8 0.8], ...
  'Style','text', ...
  'String','   seconds remaining until analysis continues', ...
  'Tag','headline1');

hui1 = uicontrol('Parent',fh_bscore, ...
  'Units','normalized', ...
  'HorizontalAlignment','left', ...
  'Position',[.5 .2 .18 .13], ...
  'FontSize',18, ...
  'Fontweight','bold',...
  'style','pushbutton',...
  'String','OK & continue', ...
  'TooltipString','press to proceed before countdown ends',...
  'ListboxTop',0, ...
  'Tag','goBttn', ...
  'callback', 'delete(gcbf)');
hui2 = uicontrol('Parent',fh_bscore, ...
  'Units','normalized', ...
  'HorizontalAlignment','left', ...
  'Position',[.5 .06 .18 .13], ...
  'FontSize',18, ...
  'Fontweight','bold',...
  'style','pushbutton',...
  'String','not OK, break', ...
  'TooltipString','press to halt analysis',...
  'ListboxTop',0, ...
  'Tag','nogoBttn', ...
  'callback', 'evalin(''caller'',''brk=1'');');

% give 4 s time to OK or veto
% in the next version try uiwait & uiresume
i=8;
brk=0;
while i>0 && ishandle(fh_bscore) && ~brk
  if ishandle(fh_bscore),
    set(hui0,'string',[int2str(i) ' seconds remaining until analysis continues']);
    pause(1);
  end
  i=i-1;
end
if brk, 
  error('user break requested');  
end
if ishandle(fh_bscore), 
  delete(fh_bscore); 
end

% *** set time frame:
etsl(:,etslc.tsCol)=etsl(:,etslc.tsCol)+b_delta;



% ----- LOCAL FUNCTIONS --------- LOCAL FUNCTIONS --------- LOCAL FUNCTIONS ----
% ----- LOCAL FUNCTIONS --------- LOCAL FUNCTIONS --------- LOCAL FUNCTIONS ----


function etsl=sspulse2tsl(d,trigThresh,si,mfac,tag,trigTauPts)
% sspulse2tsl = 'Sampled Series of Pulses 2 Time Stamp List' generates
% extended time stamp list (etsl) from single channel data D containing 
% trigger pulses. Each trigger event will be marked according to its type
% (= trigger level) in a column of etsl. Input var 'tag' is used for this.
% Both D and the thresholds will be multiplied by input var mfac. This way,
% the task of negative-going pulses and negative thresholds can be dealt
% with without change of code. !NOTE: mfac*trigThresh MUST BE SORTED IN
% ASCENDING ORDER!
global DS
etslconst;

if mfac ~=1.0
  d=d*mfac;
  trigThresh=trigThresh*mfac;
end
if any(diff(trigThresh)<=0)
  error('trigger thresholds not properly sorted');
end
[n1 n2]=size(d);
nTrigLevel=length(trigThresh);
if length(tag)~=nTrigLevel, error('sspulse2tsl: check input vars'); end
ind=cell(1,nTrigLevel);
% start with lowest positive threshold, and work way up
for i=1:nTrigLevel
  if i==1
    % 1st run
    ind{i}=find(d>=trigThresh(i,1));
  else
    ind{i}=d(ind{i-1})>=trigThresh(i,1);
    ind{i}=ind{i-1}(ind{i});
  end
end
if isempty(ind{1})
  etsl=[];
else
  % round2: 
  for i=nTrigLevel:-1:2
    if ~isempty(ind{i})
      % extend the length of the trigger level found by trigTauPts on either side of pulse
      % ind MUST be column vectors
      if size(ind{i},2)>1,
        warning('internal problem:''ind'' must be a column vector - adjusting now');
        ind{i}=makecol(ind{i});
      end
      tmpi=find(diff(ind{i},1,1)>1);
      % the stop points of trigger plateaus - include last point explicitly and extend by trigTau
      stopPt=ind{i}([tmpi; end])+trigTauPts;
      % the start points of trigger plateaus - include last point explicitly and extend by trigTau
      startPt=ind{i}([1; tmpi+1])-trigTauPts;
      ind{i}=[];
      for ti=1:size(startPt,1)
        ind{i}=cat(1,ind{i},(startPt(ti):stopPt(ti))');
      end
      % from ALL trigger levels below the current one delete indices belonging to the
      % current one (including flanks)
      for j=i-1:-1:1
        ind{j}=setdiff(ind{j},ind{i});
      end
    end
  end
  % next step: generate tsl
  for i=1:nTrigLevel
    if isempty(ind{i})
      tsl{i}=[];
    else
      % get rid of points representing the plateau of trigger pulse, keeping only the first
      % (dead time: 2 sampling intervals). On the occasion, 
      % a) correct offset (trigTauPts) imposed above (BUT NOT FOR LOWEST LEVEL)
      % b) convert time stamps to ms
      tsl{i}=si*1e-3*(tsldeadt(ind{i},2,'include_1st',1)+(i>1)*trigTauPts);
    end
  end
  % finally, generate extended tsl from all tsl
  etsl=[];
  for i=1:nTrigLevel
    if ~isempty(tsl{i})
      tmp=tsl{i};
      tmp(:,etslc.tagCol)=tag(i);
      etsl=cat(1,etsl,tmp);
    end
  end
  % do not sort etsl because they will be merged in main function anyway
end

