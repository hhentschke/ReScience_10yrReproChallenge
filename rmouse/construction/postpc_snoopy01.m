function helper01
% ** function helper01
% an open script probing this and that:
% reads small excerpts of streams corresponding to the time stamps of specific populations
% identified in principal components analysis and plots some of their characteristics

global DS AP WP R

rmouse_ini;
rawCh=rmouse_chan;

strmType={'delta','theta','thetaEnv','thetaLo','thetaLoEnv','thetaHi','thetaHiEnv','gamma','gammaEnv'};
strmType={'delta','theta','gamma','gammaEnv'};
nStrms=length(strmType);

if isempty(strfind(DS.dpath,':')), DS.dpath=[WP.rootPath DS.dpath]; end
if isempty(strfind(AP.strmDir,':')), AP.strmDir=[WP.rootPath AP.strmDir]; end

% get si
abfi=abfinfo([DS.dpath '\' DS.abfFn '.abf'],'verbose',0);      
si=abfi.si;

% determine the excerpt to read from stream file (in points): 
% - the first and last intervals' midpoints in seconds from beginning of recording
intv=[R.tb(min(cat(1,WP.pop(:).ix)),1)  R.tb(max(cat(1,WP.pop(:).ix)),1)];
% - convert intv to points and extend by half of AP.ppseg on either side
halfWinPts=(AP.ppSeg/2);
if halfWinPts~=fix(halfWinPts), 
  error('AP.ppSeg is uneven');
end
intvPts=cont2discrete(intv,si*1e-6,'intv',1)+[-halfWinPts+1 halfWinPts+1];

% now compute, for each population, the observations' indices into the excerpt, taking
% into account that we have an offset due to the fact that time zero for R.tb and the
% excerpt are different
for g=1:WP.nPop
  % the start points
  tmp=cont2discrete((R.tb(WP.pop(g).ix,1))',si*1e-6,'intv',0)-(halfWinPts-1)-intvPts(1)+1;
  % expand columns so we have a column array with ready-to-use indices
  WP.pop(g).intvIx=repmat(tmp,AP.ppSeg,1)+repmat((0:AP.ppSeg-1)',1,length(tmp));
end
  

% generate one var for each stream, this is most flexible for all sorts of plots,
% including pllplot overlays
for ci=length(AP.LFPccInd):-1:1
  for six=1:length(strmType)
    % sic: rawCh(ci)
    eval(['tmpD=strmread([AP.strmDir ''\'' rawCh(ci).' strmType{six} 'Fn],''intv'',intvPts,''verbose'',0);' ]);
    for g=1:WP.nPop
      tmpr=evdeal(tmpD(WP.pop(g).intvIx),'idx',{'rms','area'});
      eval(['WP.pop(g).' strmType{six} 'rms(1:WP.pop(g).nObs,ci)=tmpr.rms'';']);
    end
  end
end

popmean=zeros(nLFPCh,WP.nPop);
for six=1:length(strmType)
  figure
  for g=1:WP.nPop
    eval(['popmean(:,g)=mean(WP.pop(g).' strmType{six} 'rms,1)'';']);
    leg{g}=['pop ' int2str(g)];
  end
  plot(popmean,'o-');
  legend(leg);
  title(strmType{six});
end

    

% counterpart to strmwrite
function d=strmread(varargin)
global DS
d=i16load(varargin{:})/2^15*DS.nsRng;
