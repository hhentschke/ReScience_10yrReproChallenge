function movie_ccShortInterval01
% demonstrates the rapid dynamics of theta cc lag changes
global DS AP 

% subplots: raw data (all ch), raw cc, (contour plot cc), matrix CC lag, matrix CC peak,
% time course cc lag

% improvements:
% - set(<plot handle>,'xdata',... instead of repeated plot commands
% - y axis labels in terms of periods
% - choice of channels, dist of recording sites, etc. automatic/more user-friendly
% - implement 'princ' phase lags

% -----------------------------------------------------------------------------
%                          I. PRELIMINARIES
% -----------------------------------------------------------------------------
plotMd=1;
rmouse_ini;
mouse=1;
switch mouse
  case 1
    % --- session-specific settings
    intv=[22.5 26]; % a nice phase jump in there!
    intv=[190 194]; % a nice phase jump in there!    
    intv=[202 207]; % fat phase jump
    intv=[148 153]; % corresponds to movie sequence where mouse rears up
    intv=[1000 1005]; % a nice phase jump in there!

    intv=[120 123]; % very interesting excerpt with one nice sequence of one generator overtaking the other
    intv=[126 131]; % corresponds to movie sequence where mouse finally decides to get up and walk around
    
    cd([WP.rootPath '\beta3_wtko\wt0001_04708']);
    a001_r1;
    % ** here is the place to reduce the subset of channels to be used in this
    % routine **
    AP.rawChAnNm=AP.rawChAnNm(3:end-2);
    % lag(s) of peak within CC which to pick, expressed in ms
    peakLag=[20];
    peakLag=[40];    
    ccThresh=.5,
    
  case 2
    intv=[300 302];
    cd([WP.rootPath '\beta3_wtko\ko2445_03214']);
    r1;
    % ** here is the place to reduce the subset of channels to be used in this
    % routine **
    AP.rawChAnNm=AP.rawChAnNm(3:end-2);
    AP=rmfield(AP,'OPT_thgaeCCAttractLag');
    peakLag=[70];    
    ccThresh=.5,    
  
  case 3
    intv=[1500 1505]+60;
    cd([WP.rootPath '\beta3_wtko\ko2483_03303']);
    r1;
    % ** here is the place to reduce the subset of channels to be used in this
    % routine **
    AP.rawChAnNm=AP.rawChAnNm(1:end-5);
    peakLag=[60];    
    ccThresh=.5,    
  otherwise
    error('sdfsdf');
end

% --- channels
rmouse_APcheck;
rawCh=rmouse_chan;
if isempty(strfind(DS.dpath,':')), DS.dpath=[WP.rootPath DS.dpath]; end
if isempty(strfind(AP.strmDir,':')), AP.strmDir=[WP.rootPath AP.strmDir]; end

% index (into AP.rawChAnNm) for the channel pair from which to show raw CC results:
% channel xx um more dorsal of principal (=lm) channel
dDist=.5;

[doff,ccRefChInd]=min(abs(WP.elx(AP.LFPccInd)+dDist));
ccRefChIdx=AP.LFPIdx(ccRefChInd);
dDist=diff(WP.elx([ccRefChIdx AP.pcIdx]));
subChInd=[ccRefChInd AP.pcInd];

% --- stream type (one only) and appropriate windows
strmType={'theta'};


nStrms=length(strmType);
AP.ppSeg=                     256;
AP.dftOlapPts=                231;
AP.ccLagPts=                    150;   
AP.ccScaleOpt=                 'coeff';'coeff_ub';
AP.thccw=                      [-145 145];
AP.gaccw=                      [-40 40];
AP.cca=                        [-144 144];


% --- load data
% generate one var for each stream, this is most flexible for all sorts of plots,
% including pllplot overlays
for ci=nCh:-1:1
  % start by loading raw data from abf file, thus obtaining si
  if ci==nCh
    [rawD,si]=abfload([DS.dpath '\' DS.abfFn '.abf'],'start',intv(1),'stop',intv(2),'channels',AP.rawChAnNm);        
    if DS.rawSignalInverted
      rawD=-1*rawD;
    end
    intvPts=cont2discrete(intv*1e6,si,'intv',1);
  end
  for six=1:length(strmType)
    eval(['d(:,ci)=strmread([AP.strmDir ''\'' rawCh(ci).' strmType{six} 'Fn],''intv'',intvPts,''verbose'',0);' ]);
%     env=strmread([AP.strmDir '\' rawCh(ci).thetaHiEnvFn],'intv',intvPts,'verbose',0);
%     d(:,ci)=d(:,ci)./env;
  end
  ndPts=size(d,1);
end

% --- subdivide data in intervals
[subInt,subInt_pts]=mkintrvls([1 ndPts],'resol',1,'ilen',AP.ppSeg,'olap',AP.dftOlapPts,'border','skip');
nSubInt=size(subInt,1);

% --- graphics
% standard settings for raw plot
labelscale('fontSz',10,'scaleFac',1,'lineW',.25,'markSz',4); 
% figdir='d:\projects\madison\rmouse\sfn2004\ps\';
% printas='-dpsc2';[];
ftag='ccEvol';
figha=findobj('tag',ftag);
if isempty(figha), figha=figure;
else  figure(figha);
end
tmpScrSz=get(0,'Screensize')*.8+.1*rand;
tmpScrSz([1 2])=tmpScrSz([1 2])+tmpScrSz([3 4])*.10;
set(figha,'position',tmpScrSz,'tag',ftag,'name',ftag,...
  'color',[.9 .9 1],'numbertitle','off');
orient landscape;
clf;

% -----------------------------------------------------------------------------
%                          II. INITIAL COMPUTATION & PLOTS
% -----------------------------------------------------------------------------
xmargin=.05;
ymargin=.05;
% lags of cc in ms
lag_ms=discrete2cont([-AP.ccLagPts:AP.ccLagPts]+1,si*.001);

% --- raw data - 2/3 of figure
sph(1)=subplot('position',[xmargin 1/3+ymargin 2/3-2*xmargin 2/3-2*ymargin]);
hold on
[yl,dy]=pllplot(d,'si',si,'noscb',1);
% first plot is always =last child
c=get(gca,'children');
% make all gray
set(c(end-nCh+1:end),'color',[.5 .5 .5]);

% --- time course of lag of CC pair right below
sph(2)=subplot('position',[xmargin ymargin 2/3-2*xmargin 1/3-2*ymargin]);
axis manual, hold on;
set(gca,'xlim',subInt([1 end],1)+diff(subInt(1,:))*.5,'ylim',[-10 80]);
box on
xlabel('time (ms)');
ylabel('lag (ms)');
% upper right: peak lag and amplitude of CC(ref,x)
%...



% lower right: current CC waveform
sph(3)=subplot('position',[2/3+xmargin ymargin 1/3-2*xmargin 1/3-2*ymargin]);
axis manual, hold on;
set(gca,'xlim',AP.thccw,'ylim',[-1 1]);
xlabel('lag (ms)');
ylabel('crosscorr');

% -----------------------------------------------------------------------------
%                          II. LOOP COMPUTATIONS & PLOTS
% -----------------------------------------------------------------------------

% current plot version of d
cpd=d+nan;
cpd([1 end],:)=d([1 end],:);
% crosscorr waveforms
CC=zeros(2*AP.ccLagPts+1,nSubInt);
% CC lag
ccPeakLag=zeros(nSubInt,1);
% instantaneous theta freq of both channels (as estimated by number of peaks per time)
iThetaFreq=zeros(nSubInt,2);
% where to fish for peaks
tP=cont2discrete(peakLag,si*.001,'intv',0)+AP.ccLagPts;

for subi=1:nSubInt
  % --- all calculations:
  % - cc  
  [cc,lag]=xxcorr(d(subInt(subi,1):subInt(subi,2),subChInd(1)),...
    d(subInt(subi,1):subInt(subi,2),subChInd(2)),...
    AP.ccLagPts,AP.ccScaleOpt);
  CC(:,subi)=cc;
  % - lag of cc peak
  tmpr=evdeal(cc,'idx',{'closepeak'},'tP',tP);
  ccPeakLag(subi)=lag_ms(tmpr.tlPosPeakT);
  tlPosPeak=tmpr.tlPosPeak;

%   % - instantaneous theta freq 
%   % find max peak; this is the start lag
%   tmpr=evdeal(d(subInt(subi,1):subInt(subi,2),subChInd),'idx',{'minmaxpeak'});  
%   % sine with variable freq 
%   ft_ = fittype('a1*cos(a2*x+a3)' ,...
%     'dependent',{'y'},'independent',{'x'},...
%     'coefficients',{'a1','a2','a3'});
%   fo_ = fitoptions('method','NonlinearLeastSquares');
%   for gg=1:2
%     st_ = [.5 8.5   tmpr.maxPeakT(gg)/AP.ppSeg*discrete2cont(AP.ppSeg,si*1e-6)'*2*pi];
%     set(fo_,'Startpoint',st_);
%     fifi=fit(discrete2cont(1:AP.ppSeg,si*1e-6)'*2*pi,detrend(d(subInt(subi,1):subInt(subi,2),subChInd(gg))),ft_ ,fo_);    
%     tmpv=coeffvalues(fifi);
%     iThetaFreq(subi,gg)=tmpv(2);
%   end

  
  if plotMd
    % --- raw plot of excerpt
    cpd(subInt(subi,1):subInt(subi,2),:)=d(subInt(subi,1):subInt(subi,2),:);
    subplot(sph(1));
    pllplot(cpd,'si',si,'spacing','fixed','dy',dy,'ylim',yl,'noscb',1);
    c=get(gca,'children');
    nChildren=length(c);
    % identify the pointers to the two CC channel segments:
    % - anything added to a plot is added to the HEAD of the array of children, shifting 
    %   the position of all other items to the right
    % - pllplot creates one child per channel (excerpt)
    % -> both facts above mean that the order of channels in c is reversed with respect to
    %   d
    % - nothing except the traces is plotted here, so the number of added children is = nCh
    % So, in essence we only need to reverse subChInd
    subChRawPlotInd=nCh-subChInd+1;
    % all in black..
    set(c(1:nCh),'color',[0 0 0],'linewidth',2);
    % except the pair from which cc is calculated
    set(c(subChRawPlotInd),'color',[.8 .1 .2],'linewidth',2);
    % overwrite data just plotted with nans
    cpd(subInt(subi,1):subInt(subi,2),:)=nan;  
    % --- waveform
    subplot(sph(3));
    ccph=plot(lag_ms,cc,'-');
    set(ccph,'linewidth',2,'color','r');
    % --- lag
    % plot lags for peak crosscorrelations smaller than a certain threshold in pale color
    if tlPosPeak>=ccThresh
      pCol='r';
    else
      pCol=[1 .7 .7];
    end
    subplot(sph(2));
    ph=plot(mean(subInt(subi,:)),ccPeakLag(subi),'o');
    set(ph,'color',pCol);

%     plot(mean(subInt(subi,:)),iThetaFreq(subi,1),'g-o');
%     plot(mean(subInt(subi,:)),iThetaFreq(subi,2),'m-s');    
    
    
    drawnow; 
    % M(subi)=getframe(gcf);
    
    % pause(.1);
    % things to do right after pause:
    % - delete raw excerpt plot
    delete(c(1:nCh));  
    % - diminish cc
    set(ccph,'linewidth',.5,'color',[.5 .5 .5]);
  end  
end

disp('done');
% movie(M);

% -----------------------------------------------------------------------------
%                          LOCAL FUNCS
% -----------------------------------------------------------------------------

% counterpart to strmwrite
function d=strmread(varargin)
global DS
d=i16load(varargin{:})/2^15*DS.nsRng;

