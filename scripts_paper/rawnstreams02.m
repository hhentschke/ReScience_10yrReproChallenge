function rawnstreams02
% generates plots of raw data and selected streams
% raw & theta overlaid, gamma +gammaEnv overlaid; one channel
global DS AP WP

ch={'IN 11'};
intv=[999.5 1000.1];
intv=[99.6 100.4]-.45;
intv=[524.1  524.9]; % almost perfect similarity
intv=[524.1  525.5]; % almost perfect similarity, then disruption
intv=[99.5 100.6]+44; % a realistic picture: similarity, but not perfect
intv=[99.55 100.8]+45; % a realistic picture: similarity, but not perfect

intv=[99.45 100.55]+45; % a realistic picture: similarity, but not perfect

rmouse_ini;

cd([WP.rootPath '\beta3_wtko\wt0001_04708']);
a001_r1;

labelscale('fontSz',8,'scaleFac',.7,'lineW',.5,'markSz',8); 
ornt='portrait';
figdir='d:\projects\hippocampus\rmouse\paper_atropine\rawFig\';
% figdir='';
printas='-dpsc2';
printas=[];

figName=mfilename;

nCh=length(ch);
% local indices to channels
chInd=[];
for i=1:nCh
  chInd=[chInd strmatch(ch{i},AP.rawChAnNm,'exact')];
end
if isempty(chInd), error('check channel names'); end

strmType={'theta','gamma','gammaEnv'};
nStrms=length(strmType);

for i=1:length(AP.rawChAnNm)
  rawCh(i).nm=AP.rawChAnNm{i};
  tmpNm=[DS.abfFn '_' AP.rawChAnNm{i}];
  rawCh(i).thetaFn=[tmpNm '_theta.i16'];
  rawCh(i).thetaEnvFn=[tmpNm '_thetaEnv.i16'];
  rawCh(i).thetaLoFn=[tmpNm '_thetaLo.i16'];
  rawCh(i).thetaLoEnvFn=[tmpNm  '_thetaLoEnv.i16'];                  
  rawCh(i).thetaHiFn=[tmpNm '_thetaHi.i16'];
  rawCh(i).thetaHiEnvFn=[tmpNm  '_thetaHiEnv.i16'];                  
  rawCh(i).gammaFn=[tmpNm  '_gamma.i16'];        
  rawCh(i).gammaEnvFn=[tmpNm  '_gammaEnv.i16'];                
  rawCh(i).rippleFn=[tmpNm  '_ripple.i16'];        
  rawCh(i).rippleEnvFn=[tmpNm  '_rippleEnv.i16'];                
  rawCh(i).deltaFn=[tmpNm '_delta.i16'];  
end
clear tmpNm

if isempty(strfind(DS.dpath,':')), DS.dpath=[WP.rootPath DS.dpath]; end
if isempty(strfind(AP.strmDir,':')), AP.strmDir=[WP.rootPath AP.strmDir]; end

% generate one var for each stream, this is most flexible for all sorts of plots,
% including pllplot overlays
for ci=length(chInd):-1:1
  % start by loading abf file, thus obtaining si
  if ci==length(chInd)
    [rawD,si]=abfload([DS.dpath '\' DS.abfFn '.abf'],'start',intv(1),'stop',intv(2),'channels',ch);        
    if DS.rawSignalInverted
      rawD=-1*rawD;
    end
    intvPts=cont2discrete(intv*1e6,si,'intv',1);
    if diff(intvPts)+1>size(rawD,1)
      intvPts(end)=intvPts(end)-1;
    end
  end
  for six=1:length(strmType)
    eval([strmType{six} 'D(:,ci)=strmread([AP.strmDir ''\'' rawCh(chInd(ci)).' strmType{six} 'Fn],''intv'',intvPts,''verbose'',0);' ]);
  end
end

% thetaD=bafi(rawD,si,[8 12],'rs',50);
% gammaD=bafi(rawD,si,[15 90],'rs',50);
% gammaEnvD=abs(hilbert(gammaD));

figure(1), clf, orient(ornt)

dD=permute(cat(3,rawD,gammaD),[1 3 2]);
[n1 n2 n3]=size(dD);
dD=reshape(dD,[n1 n3*n2 1]);
[yl,dy]=pllplot(dD,'si',si,'spacing','maxmin','noscb',1);
hold on
% first plot==last child
c=get(gca,'children');
set(c([1 2]),'color',[.4 .4 .4]);

dD=permute(cat(3,thetaD,gammaEnvD),[1 3 2]);
[n1 n2 n3]=size(dD);
dD=reshape(dD,[n1 n3*n2 1]);
% scalebar
[yl,dy]=pllplot(dD,'si',si,'spacing','fixed','dy',dy,'ylim',yl,'noscb',0);
% different colors for different streams
c=get(gca,'children');
% first plot==last child
set(c([4 5]),'linewidth',get(c(1),'linewidth')*3);
% set(c([1 3 5]),'color',[.25 .7 .25]);
% set(c([7 9 11]),'color',[.75 .75 .75]);

rexy('ax',gca,'xfac',.8,'yfac',.60);

if ~isempty(printas), 
  print(printas,[figdir figName]); 
end


% counterpart to strmwrite
function d=strmread(varargin)
global DS
d=i16load(varargin{:})/2^15*DS.nsRng;
