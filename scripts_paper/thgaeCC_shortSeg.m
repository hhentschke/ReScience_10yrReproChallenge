function thgaeCC_shortSeg
% generates 
% - plot of theta, gamma +gammaEnv overlaid for figure depicting short-term
% dynamics of theta gammaEnv cc
% - contour plot theta gammaEnv cc
% - cc func theta gammaEnv cc

global DS AP 

if 1
  % wt0001
  ch={'IN 8'};

  intv=[0 3]+4; % nice large disturbance in there
  markSegInd=[10 20];

  intv=[0 3]+16.98; % OK example
  markSegInd=[5 10];

  intv=15*60+[0 3]+37.2; % OK example
  markSegInd=[7 15];
  
  intv=[171.3 172.4];
  markSegInd=[5 10];

else
  % wt0002
  ch={'IN 10'};
  ch={'IN 8'};
  
  intv=[0 3]+59.98; % so so
  markSegInd=[12 15];
  
  intv=14*60+[0 3]+13.5; % good
  markSegInd=[18 21];
  
  intv=14*60+[0 2.5]+13.9; % good
  markSegInd=[14 17];
end

% --- settings figure 1
labelscale('fontSz',8,'scaleFac',.75,'lineW',.25,'markSz',8); 
ornt='tall';
figdir='d:\projects\rmouse\paper_atropine\rawFig\';
% figdir='';
printas='-dpsc2';
printas=[];

figName=mfilename;



% ----- load streams ----------------------------------------------------------- 
rmouse_ini;

cd([WP.rootPath '\beta3_wtko\wt0001_04708']);
a001_r1;

% cd([WP.rootPath '\beta3_wtko\wt0002_04707']);
% a000_r1;

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

% ----- peliminaries for cc  --------------------------------------------
AP.ppSeg=                     128;
AP.dftOlapPts=                 28;
AP.ccLagPts=                    122;   

% intervals to use for computation of cc
[intrvls,intrvls_pts,nilen,nilen_pts]=mkintrvls(intvPts-intvPts(1)+1,...
  'resol',1,'ilen',AP.ppSeg,'olap',AP.dftOlapPts,'border','skip');
nInt=size(intrvls,1);


% --- plot streams -------------------------------------------------------------- 
figure(1), clf, orient(ornt)
subplot(4,1,1)

dD=permute(cat(3,thetaD,gammaD),[1 3 2]);
[n1 n2 n3]=size(dD);
dD=reshape(dD,[n1 n3*n2 1]);
[yl,dy]=pllplot(dD,'si',si,'spacing','maxmin','noscb',1);
hold on
% first plot==last child
c=get(gca,'children');
set(c([1 2]),'color',[.6 .6 .6]);

dD=permute(cat(3,thetaD,gammaEnvD),[1 3 2]);
[n1 n2 n3]=size(dD);
dD=reshape(dD,[n1 n3*n2 1]);
% no scalebar here because we want the same time scale as in contour plot
% below
[yl,dy]=pllplot(dD,'si',si,'spacing','fixed','dy',dy,'ylim',yl,'noscb',1);
% different line widths for different streams
c=get(gca,'children');
minLineWid=get(c(1),'linewidth');
% first plot==last child
set(c([3 4]),'linewidth',minLineWid*1.5);
% set(c([4 5]),'linewidth',minLineWid*1.5);

% mark selected segments
mIx=[];
for segi=markSegInd
  mIx=union(mIx,intrvls(segi,1):intrvls(segi,2));
end
dD(setdiff(1:size(dD,1),mIx),:)=nan;
[yl,dy,nada,ph]=pllplot(dD,'si',si,'spacing','fixed','dy',dy,'ylim',yl,'noscb',1);
% c=get(gca,'children');
% set(c([4 5]),'linewidth',minLineWid*1.5);

set(ph,'linewidth',minLineWid*5);

% --- fake subplot 1 for scale bar ----------------------
subplot(4,1,2)
[yl,dy]=pllplot(dD,'si',si,'spacing','fixed','dy',dy,'ylim',yl,'noscb',0);

% --- compute cc, label segments, make contour plot ----------------------
% preallocate 
segCC=repmat(nan,2*AP.ccLagPts+1,nInt);
% cc
for segi=1:nInt
  [segCC(:,segi),lag]=xcorr(detrend(thetaD(intrvls(segi,1):intrvls(segi,2)),'constant'),...
    detrend(gammaEnvD(intrvls(segi,1):intrvls(segi,2)),'constant'),...
    AP.ccLagPts,'coeff');
end

subplot(5,1,4)
t=discrete2cont(round(mean(intrvls,2)),si*.001,'intv',0)*.001;
colormap gray
[c,h]=contourf(t,...
  lag,segCC,[-1:.1:1]);
set(h,'linestyle','none');
set(gca,'xlim',[0 t(end)+diff(t(1:2))*.5]);
line(get(gca,'xlim'),[0 0],'linestyle','--','color',[.6 .6 .6]);

% % --- fake subplot 2 for color bar ----------------------
% subplot(5,1,5)
% [c,h]=contourf(t,...
%   lag,segCC,[-1:.1:1]);
% colorbar('EastOutside');


if ~isempty(printas), 
  print(printas,[figdir figName]); 
end


% --- figure 2 for cc funcs ----------------------
% --- figure 2 for cc funcs ----------------------
figure(2), clf, orient portrait
labelscale('fontSz',8,'scaleFac',.35,'lineW',.5,'markSz',8);
subplot(2,2,1)
% all
ph=plot(lag,segCC,'k-');
set(ph,'color',[.6 .6 .6]);
hold on
lw=get(ph(1),'linewidth');
% selected
ph=plot(lag,segCC(:,markSegInd),'k-');
set(ph,'linewidth',lw*2);
axis tight
set(gca,'ylim',[-1.05 1.05]);
xl=get(gca,'xlim');
% lines
lh=line(xl,[0 0],'linestyle','--','linewidth',lw*.5,'color','k');
lh=line([0 0],get(gca,'ylim'),'linestyle','--','linewidth',lw*.5,'color','k');
set(gca,'xtick',[fliplr([0:-100:xl(1)]) [100:100:xl(2)]]);
set(gca,'ytick',[-1:.5:1]);
xlabel('Lag (ms)');

if ~isempty(printas), 
  print(printas,[figdir figName '_2']); 
end


% --- local func --------------------------------------------------------------
% --- local func --------------------------------------------------------------

% counterpart to strmwrite
function d=strmread(varargin)
global DS
d=i16load(varargin{:})/2^15*DS.nsRng;
