function rawnstreams01
% generates plots of raw data and selected streams
% result section, comparison ctrl - atropine
global DS AP 

labelscale('fontSz',8,'scaleFac',.8,'lineW',.25,'markSz',8); 
ornt='landscape';
ornt='tall';

figdir='d:\projects\rmouse\paper_atropine\rawFig\';
% figdir='';
printas=[];
% printas='-dpsc2';
figName=mfilename;
figure(1), clf, orient(ornt)

% standard for to-be-used plots
exc=[-1.25 1.25];
exc=[-1.5 1.5];
% good for snooping
% exc=[-4 4];


rmouse_ini;

exprmnt={'\beta3_wtko\wt2654_03304','\beta3_wtko\wt2654_03304'};
exprmnt={'\beta3_wtko\wt0001_04708','\beta3_wtko\wt0001_04708'};
exprmnt={'\beta3_wtko\wt0003_04730','\beta3_wtko\wt0003_04730'};


% loop over experiments 
for k=1:length(exprmnt)
  AP=[]; DS=[];
  switch exprmnt{1}
    case '\beta3_wtko\wt0001_04708'
      switch k
        case 1
          cd([WP.rootPath exprmnt{k}]);
          a001_r1;
          ch=AP.rawChAnNm(6:12);
          intv=599.4+exc;
        case 2
          cd([WP.rootPath exprmnt{k}]);
          a003_r1;
          ch=AP.rawChAnNm(6:12);
          intv=992.2+exc; % soso
          intv=1010.25+exc; % not bad
          intv=1485.2+exc;      % good

          intv=311.2+exc;      % good
          intv=536.4+exc;      % good

      end % switch:k
    case '\beta3_wtko\wt0003_04730'
      switch k
        case 1
          cd([WP.rootPath exprmnt{k}]);
%           a002_r1;
%           ch=AP.rawChAnNm(3:9);
%           intv=29.25+exc; % not convincing
%           intv=745.7+1.5+exc; % better 
%           intv=1588.4+exc; % better 

          a003_r1;
          ch=AP.rawChAnNm(3:9);
          intv=365.1+exc; % 
          
          
        case 2
          cd([WP.rootPath exprmnt{k}]);
          a005_r1;
          ch=AP.rawChAnNm(3:9);
          intv=430.0+exc; % soso
          intv=1564.2+exc; % soso

%           a004_r1;
%           ch=AP.rawChAnNm(3:9);
%           intv=405+exc; % soso          
%           intv=571+exc; % soso          
          
%           0.0093
%           0.1562
%           0.3135
%           0.3356
%           0.3619
%           0.3912
%           0.4059
%           0.4610
%           0.4770
%           0.5279
%           0.5552
%           0.5736
%           0.9965
%           1.5181
      end % switch:k
    case '\beta3_wtko\wt2654_03304'
      switch k
        case 1
          cd([WP.rootPath exprmnt{k}]);
          a000_r1;
          ch=AP.rawChAnNm(2:8);
          intv=1901+exc; % not convincing
          intv=1790+exc; % random!
        case 2
          cd([WP.rootPath exprmnt{k}]);
          a001_r1;
          ch=AP.rawChAnNm(2:8);
          intv=1990+exc; % random!
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
    
    subplot(3,1,k)
    if k==1
      [yl,dy]=pllplot(rawD,'si',si,'noscb',1,'noplot',1,'spacing','fixed','dy',1.3);
    end
    pllplot(rawD,'si',si,'spacing','fixed','dy',dy,'ylim',yl);
    rexy('ax',gca,'xfac',.75,'yfac',1);

end

if ~isempty(printas), 
  print(printas,[figdir figName]); 
end


% ---------------------------------------------------------------------

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

