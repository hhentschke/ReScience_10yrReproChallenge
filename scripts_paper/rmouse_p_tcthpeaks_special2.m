% ------- plots of time course:
intv=AP.rawExcerpt;
boe=intv(1)*1000;
eoe=intv(2)*1000;
ch=AP.rawChAnNm;
% local indices to channels
chInd=[];
for i=1:nCh
  chInd=[chInd strmatch(ch{i},AP.rawChAnNm,'exact')];
end

strmType={'gammaEnv'};

% generate one var for each stream, this is most flexible for all sorts of plots,
% including pllplot overlays
for ci=length(chInd):-1:1
  % start by loading abf file, thus obtaining si
  if ci==length(chInd)
    if exist([DS.dpath '\' DS.abfFn '.mat'],'file')
      [rawD,si]=matDload([DS.dpath '\' DS.abfFn '.mat'],'start',intv(1),'stop',intv(2),'channels',ch);
    else
      [rawD,si]=abfload([DS.dpath '\' DS.abfFn '.abf'],'start',intv(1),'stop',intv(2),'channels',ch);        
    end
    if DS.rawSignalInverted
      rawD=-1*rawD;
    end
    intvPts=cont2discrete(intv*1e6,si,'intv',1);
  end
  for stIx=1:length(strmType)
    eval(['fNm=[AP.strmDir ''\'' rawCh(chInd(ci)).' strmType{stIx} 'Fn];']);
    if exist(fNm,'file')
      eval([strmType{stIx} 'D(:,ci)=strmread([AP.strmDir ''\'' rawCh(chInd(ci)).' strmType{stIx} 'Fn],''intv'',intvPts,''verbose'',0);' ]);
    end
  end
end

figure(1), clf
colormap bone
subplot(4,1,1)
pickFac=5;
gammaEnvD=lofi(gammaEnvD,si,20,'rs',30,'pickf',pickFac);
si=si*5;
% plot vs time in s
[c,cph]=contourf([1:size(gammaEnvD,1)]/1e6*si,1:7,gammaEnvD',20);
set(cph,'linestyle','none');
set(gca,'ydir','rev')
cl=get(gca,'clim');
set(gca,'clim',[cl(1) cl(2)*.9]);
axis off
colorbar

% ====================================================================

labelscale('fontSz',8,'scaleFac',1.0,'lineW',.25,'markSz',4); 

% -- preps

markSz=11;
mCol='y';

pcIdx=AP.pcIdx;
pcInd=AP.pcInd;
% -- calculations
load(WP.thetaPFn);
pp_nppc=repmat(nan,1,nLFPCh);
np_nppc=repmat(nan,1,nLFPCh);
for chInd=AP.LFPInd
  % -- positive-going peaks:
  % times of occurrence in s
  pt=.001*posPT{AP.LFPccInd(chInd)};
  tmpi=find(pt>=boe*.001 & pt<eoe*.001);
  % retain number of peaks per channel (see below)
  pp_nppc(chInd)=length(tmpi);
  pp_tpm(1:pp_nppc(chInd),chInd)=pt(tmpi);
  % amplitudes
  pp_apm(1:pp_nppc(chInd),chInd)=posPA{AP.LFPccInd(chInd)}(tmpi);
  % -- negative-going peaks:
  % times of occurrence
  pt=.001*negPT{AP.LFPccInd(chInd)};
  tmpi=find(pt>=boe*.001 & pt<eoe*.001);
  % retain number of peaks per channel (see below)
  np_nppc(chInd)=length(tmpi);
  np_tpm(1:np_nppc(chInd),chInd)=pt(tmpi);
  % amplitudes
  np_apm(1:np_nppc(chInd),chInd)=negPA{AP.LFPccInd(chInd)}(tmpi);
end
clear pt tmpi;
% pp_tpm and pp_tpm accomodate column vectors which originally had different 
% lengths; the zeroes of their zero-padding should be set to nan:
if length(isfinite(unique(pp_nppc)))~=1,
  pp_mnppc=size(pp_tpm,1);
  for chInd=AP.LFPInd
    if pp_nppc(chInd)<pp_mnppc,
      pp_tpm(pp_nppc(chInd)+1:end,chInd)=nan;
    end
  end
end
if length(isfinite(unique(np_nppc)))~=1,
  np_mnppc=size(np_tpm,1);
  for chInd=AP.LFPInd
    if np_nppc(chInd)<np_mnppc,
      np_tpm(np_nppc(chInd)+1:end,chInd)=nan;
    end
  end
end

pp_tpm=pp_tpm-intv(1);
np_tpm=np_tpm-intv(1);


% (what a heap of code for such a simple task..  )
% invert neg-going peak amplitudes
np_apm=np_apm*-1;
% negative amplitude of pos-going peaks and pos amplitude of neg-going peak 
% will be assigned smallest size: 
pp_apm(pp_apm<=0)=eps;
np_apm(np_apm<=0)=eps;
% scale each channel individually: disregard (in setting the upper limit) 
% amplitudes representing upper .5 % of points since these very likely 
% represent artifacts (which spoil scaling)

pp_co=cumh(pp_apm,.005,'p',[.995]);
np_co=cumh(np_apm,.005,'p',[.995]);

pp_co=ones(size(pp_co))*.5;
np_co=ones(size(np_co))*.5;

% replace amplitudes by marker sizes: largest amplitude of 99.5 % events=8
pp_apm=ceil(markSz*pp_apm./repmat(pp_co,size(pp_apm,1),1));
np_apm=ceil(markSz*np_apm./repmat(np_co,size(np_apm,1),1));
% no blobs from outer space, please: restrict size
pp_apm(pp_apm>markSz)=markSz;
np_apm(np_apm>markSz)=markSz;

yl=[.5 nAllLFPCh+.5];

% freeze axis, 'hold' mode, y-axis reversed so that most dorsal electrode is on top
set(gca,'NextPlot','add',...
  'XLimMode','manual','YLimMode','manual','xaxisloc','top','ydir','reverse');

for chInd=AP.LFPInd
  % --- pos peaks:
  % find all points in current interval
  tmpi=find(pp_tpm(:,chInd));
  % plot each point individually to be able to set its marker size
  for i=tmpi'
    ph=plot(pp_tpm(i,chInd),AP.LFPInd(chInd),'o');
    set(ph,'markersize',pp_apm(i,chInd),'color',mCol);
  end
  % --- neg peaks:
  % find all points in current interval
  tmpi=find(np_tpm(:,chInd));
  % plot each point individually to be able to set its marker size
  for i=tmpi'
    ph=plot(np_tpm(i,chInd),AP.LFPInd(chInd),'o');
    set(ph,'markersize',np_apm(i,chInd),'color',mCol,'markerfacecolor',mCol);
  end
end

clear pp* np* figi nspp spi ph pcol li bsbix bp_etsl

