function rawnstreams03
% generates plots of raw & theta overlaid + dots on theta peaks; one
% channel; exploring and immobile, control and atropine 

global DS AP WP

labelscale('fontSz',8,'scaleFac',.5,'lineW',.25,'markSz',3); 
ornt='portrait';

figdir='d:\projects\rmouse\paper_atropine\rawFig\';
% figdir='';
printas=[];
% printas='-dpsc2';
figName=mfilename;
figure(1), clf, orient(ornt)

% standard excerpt length
exc=[-.5 .5];
exprmnt={'\beta3_wtko\wt0001_04708','\beta3_wtko\wt0001_04708'};
yl=[-1.5 1.5];
rmouse_ini;

% loop over experiments = drug applications
for k=1:length(exprmnt)
  AP=[]; DS=[];
  switch exprmnt{1}
    case '\beta3_wtko\wt0001_04708'
      intv=[];
      switch k
        case 1
          cd([WP.rootPath exprmnt{k}]);
          a001_r1;
          ch=AP.rawChAnNm(12);
          % first row immobile, second row exploring
          intv(1,1:2)=597.7+exc;
          intv(2,1:2)=524.30+exc;          
        case 2
          cd([WP.rootPath exprmnt{k}]);
          a003_r1;
          ch=AP.rawChAnNm(12);
          intv=992.45+exc; % soso
          intv=1010.5+exc; % not bad
          intv=1485.45+exc;      % good
          intv=311.45+exc;      % good

          intv=536.65+exc;      % good
          intv(1,1:2)=535.55+exc;
          intv(2,1:2)=538.7+exc;          

      end % switch:k
  end  % switch:exprmnt

  nCh=length(ch);
  % local indices to channels
  chInd=[];
  for i=1:nCh
    chInd=[chInd strmatch(ch{i},AP.rawChAnNm,'exact')];
  end
  if isempty(chInd), error('check channel names'); end
  
  strmType={'theta'};
  nStrms=length(strmType);
  
  for i=1:length(AP.rawChAnNm)
    rawCh(i).nm=AP.rawChAnNm{i};
    tmpNm=[DS.abfFn '_' AP.rawChAnNm{i}];
    rawCh(i).thetaFn=[tmpNm '_theta.i16'];
  end
  clear tmpNm
  
  if isempty(strfind(DS.dpath,':')), DS.dpath=[WP.rootPath DS.dpath]; end
  if isempty(strfind(AP.strmDir,':')), AP.strmDir=[WP.rootPath AP.strmDir]; end
  
  for ivi=1:2
    % generate one var for each stream, this is most flexible for all sorts of plots,
    % including pllplot overlays
    for six=1:length(strmType)
      eval([strmType{six} 'D=[];']);
    end
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
        intvPts=cont2discrete(intv(ivi,:)*1e6,si,'intv',1);
      end
      for six=1:length(strmType)
        eval([strmType{six} 'D(:,ci)=strmread([AP.strmDir ''\'' rawCh(chInd(ci)).' strmType{six} 'Fn],''intv'',intvPts,''verbose'',0);' ]);
      end
    end

    rawD=bafi(rawD,si,[1 400],'rs',50);
    if size(thetaD,1)>size(rawD,1)
      thetaD=thetaD(1:size(rawD,1),:);
    end
    tmpr=evdeal(thetaD,'idx','allpeaks');

    % subplot(2,2,(k-1)*2+ivi), hold on
    subplot(2,2,(ivi-1)*2+k), hold on
    
    ph=plot((1:size(rawD,1))*si*.001,[rawD thetaD],'k');
    set(ph(1),'color',[.6 .6 .6]);
    set(ph(2),'linewidth',get(ph(1),'linewidth')*2.5);
    ph2=plot(tmpr.negPeakT{1}-1,tmpr.negPeak{1},'ko');
    set(ph2,'markerfacecolor','k');
    if k==1, axis off, end
    axis tight
    set(gca,'ylim',yl);
    % subpax(gcf)
    
    
  end
end

if ~isempty(printas), 
  print(printas,[figdir figName]); 
end


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

