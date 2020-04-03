function rawnstreams01b
% generates plots of theta and gamma streams (excerpts of what
% rawnstreams01 shows)
% all channels, ctrl and atropine
global DS AP

labelscale('fontSz',8,'scaleFac',.72,'lineW',.25,'markSz',8);
ornt='landscape';
ornt='tall';

figdir='d:\projects\hippocampus\rmouse\paper_atropine\rawFig\';
% figdir='';
printas=[];
printas='-dpsc2';
figName=mfilename;
figure(1), clf, orient(ornt)

exc=[-1.5 1.5];
exc_theta=[-.25 .25];
exc_gamma=[-.05 .05];

rmouse_ini;

exprmnt={'\beta3_wtko\wt0003_04730','\beta3_wtko\wt0003_04730'};


% loop over experiments
for k=1:length(exprmnt)
  AP=[]; DS=[];
  switch exprmnt{1}
    case '\beta3_wtko\wt0003_04730'
      switch k
        case 1
          cd([WP.rootPath exprmnt{k}]);
          a003_r1;
          ch=AP.rawChAnNm(3:9);
          intv=365.1+exc; %
          % window WITHIN this interval for theta
          intv_theta=2.7+exc_theta;
          % same for gamma
          intv_gamma=2.7+exc_gamma;

          % window WITHIN this interval for theta
          intv_theta=2.24+exc_theta;
          % same for gamma
          intv_gamma=2.36+exc_gamma;

        case 2
          cd([WP.rootPath exprmnt{k}]);
          a005_r1;
          ch=AP.rawChAnNm(3:9);
          intv=1564.2+.0+exc; % soso

      end % switch:k
  end  % switch:exprmnt
  nCh=length(ch);
  % local indices to channels
  chInd=[];
  for i=1:nCh
    chInd=[chInd strmatch(ch{i},AP.rawChAnNm,'exact')];
  end
  if isempty(chInd), error('check channel names'); end

  strmType={'theta','gamma'};
  nStrms=length(strmType);

  if isempty(strfind(DS.dpath,':')), DS.dpath=[WP.rootPath DS.dpath]; end
  if isempty(strfind(AP.strmDir,':')), AP.strmDir=[WP.rootPath AP.strmDir]; end


  ivi=1;
  for ci=length(chInd):-1:1
    % start by loading abf file, thus obtaining si
    if ci==length(chInd)
      if exist([DS.dpath '\' DS.abfFn '.mat'],'file')
        [rawD,si]=matDload([DS.dpath '\' DS.abfFn '.mat'],'start',intv(ivi,1),'stop',intv(ivi,2),'channels',ch);
      else
        [rawD,si]=abfload([DS.dpath '\' DS.abfFn '.abf'],'start',intv(ivi,1),'stop',intv(ivi,2),'channels',ch);
      end
      if DS.rawSignalInverted
        rawD=-1*rawD;
      end
    end
  end

  rawD=bafi(rawD,si,[2 400],'rs',50);
  rawD=killhum(rawD,si,60);
  rawD=killhum(rawD,si,180);

  intvPts_theta=cont2discrete(intv_theta*1e6,si,'intv',1);
  thetaD=bafi(rawD(intvPts_theta(1):intvPts_theta(2),:),si,[4 12],'rs',40);

  intvPts_gamma=cont2discrete(intv_gamma*1e6,si,'intv',1);
  gammaD=bafi(rawD(intvPts_gamma(1):intvPts_gamma(2),:),si,[40 90],'rs',40);

  % theta first
  subplot(3,2,(k-1)*2+1)
  if k==1
    [yl1,dy1]=pllplot(thetaD,'si',si,'noscb',1,'noplot',1,'spacing','fixed','dy',.55);
  end
  pllplot(thetaD,'si',si,'spacing','fixed','dy',dy1,'ylim',yl1);
  axis on
  grid on
  set(gca,'ytick',[]);
  rexy('ax',gca,'xfac',.70,'yfac',1);
  
  % gamma next
  subplot(3,2,(k-1)*2+2)
  if k==1
    [yl2,dy2]=pllplot(gammaD,'si',si,'noscb',1,'noplot',1,'spacing','fixed','dy',.23);
  end
  pllplot(gammaD,'si',si,'spacing','fixed','dy',dy2,'ylim',yl2);
  axis on
  grid on
  set(gca,'ytick',[]);
  rexy('ax',gca,'xfac',.70,'yfac',1);
  
  % for reference, raw data with marker of segments
  subplot(3,2,2*2+k), hold on
  plot(rawD(:,end),'k-'); 
  line(intvPts_theta,[0 0],'color','g','linewidth',4);
  line(intvPts_gamma,[0 0],'color','r','linewidth',4);
  axis tight


end

if ~isempty(printas),
  print(printas,[figdir figName]);
end

% ------------------------------------------------------------------
% counterpart to strmwrite
function d=strmread(varargin)
global DS
d=i16load(varargin{:})/2^15*DS.nsRng;



% loads neural data stored in matfiles (mimicks behavior of abfload)
function [d,si]=matDload(fn,varargin)
% defaults
start=0.0;
stop='e';
channels={' '};
pvpmod(varargin);
if strcmpi(channels,'a'),
  error('loading neural data from matfile: channel names must be given explicitly (probably AP.rawChAnNm is the problem)');
end
% load abfi
load(fn,'abfi');
% --- deal with start & stop times
ix(1)=cont2discrete(start*1e6,abfi.si,'intv',0);
if strcmpi(stop,'e')
  ix(2)=abfi.dataPtsPerChan;
else
  ix(2)=cont2discrete(stop*1e6,abfi.si,'intv',1);
end
% ---- load the crap
% preallocate
d=repmat(nan,diff(ix)+1,length(channels));
for chInd=1:length(channels)
  % channel names must be deblanked
  dbch=channels{chInd};
  dbch=dbch(~(int16(dbch)==32));
  load(fn,dbch);
  eval(['d(:,chInd)=double(' dbch '(ix(1):ix(2),:));']);
  eval(['clear ' dbch ';']);
end
si=abfi.si;

