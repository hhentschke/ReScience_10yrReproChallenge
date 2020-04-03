function rawnstreams04
% generates plots of raw data, comparison ctrl - iso - iso+atropine
% princ chan only
global DS AP 

labelscale('fontSz',8,'scaleFac',.8,'lineW',.25,'markSz',8); 
ornt='portrait';

figdir='d:\projects\rmouse\paper_atropine\rawFig\';
% figdir='';
printas='-dpsc2';
printas=[];
figName=mfilename;
figure(1), clf, orient(ornt)

% standard for to-be-used plots
exc=[-1.8 1.8];
% good for snooping
% exc=[-4 4];


rmouse_ini;
% list four times because we need immobile and expl for control
exprmnt={'\beta3_wtko\isoatr\wt2141_02523','\beta3_wtko\isoatr\wt2141_02523',...
  '\beta3_wtko\isoatr\wt2141_02523','\beta3_wtko\isoatr\wt2141_02523'};

cRawD=[];
% loop over experiments 
for k=1:length(exprmnt)
  AP=[]; DS=[];
  cd([WP.rootPath exprmnt{k}]);
  switch exprmnt{1}
    case '\beta3_wtko\isoatr\wt2141_02523'
      switch k
        case 1
          % control, immobile
          a000_r1;
          ch=AP.rawChAnNm(4);
          intv=803+exc; % OK
        case 2
          % control, exploring
          cd([WP.rootPath exprmnt{k}]);
          a000_r1;
          ch=AP.rawChAnNm(4);
          intv=565.0+exc; % OK
        case 3
          % iso
          a001_r1;
          ch=AP.rawChAnNm(4);
          intv=1067+exc; % 
        case 4
          % iso+atropine
          a002_r1;
          ch=AP.rawChAnNm(4);
          intv=1172+exc; % 
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
  
  for i=1:length(AP.rawChAnNm)
    rawCh(i).nm=AP.rawChAnNm{i};
    tmpNm=[DS.abfFn '_' AP.rawChAnNm{i}];
  end
  clear tmpNm
  
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
        intvPts=cont2discrete(intv(ivi,:)*1e6,si,'intv',1);
      end
    end

    rawD=bafi(rawD,si,[2 400],'rs',50);
    rawD=killhum(rawD,si,60);
    rawD=killhum(rawD,si,180);
    
    cRawD=[cRawD rawD];

end

pllplot(cRawD,'si',si,'spacing','fixed','dy',1.9);
rexy('ax',gca,'xfac',.65,'yfac',.5);


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

