function rmouse_ovrawplot(etsl,partJob)
% ** function rmouse_ovrawplot(etsl,partJob)
% generates overview plot of raw data

global DS AP WP

orient landscape
set(gcf,'DefaultLineLineWidth',.1);
remarxy('xfac',1.0,'yfac',1.05);
% philosophy: fixed number of rows and variable length of excerpt per row 
% depending on total length of data file
start=WP.boe/1000;
stop=WP.eoe/1000;
% if the number of channels is >3 pick only a subset: starting at principal 
% channel go dorsal and pick every other channel (three altogether)
if length(AP.rawChAnNm)<=3
  rawChAnNm=AP.rawChAnNm;
  local_pcIx=AP.LFPpcInd1;
else
  local_chIx=unique(round(linspace(max(1,AP.LFPpcInd1-4),AP.LFPpcInd1,3)));
  rawChAnNm=AP.rawChAnNm(local_chIx);
  local_pcIx=numel(local_chIx);
end
  
nChan=length(rawChAnNm);
% number of rows on page=number of intervals
nRow=15;
% vertical spacing of channels on plot in units of the data..
dy=2.0;
% ..and y limits, based on number of channels, their vertical spacing and 
% assumptions on signal amplitudes
yl=[-0.75-(nChan-1)  0.5]*dy;
% 4-element array to be used for patches (rectangles)
ylArr=[yl'; rot90(yl)];
% downsampling factor
sampFac=3;

etslconst;

% a local, brightened-up flavor of colors
briCol=cat(1,AP.segmentType{:,3});
% brighten colors up
briCol=max(.5,briCol);
briCol=briCol./(max(briCol,[],2)*[1 1 1])*.99;

% load
if exist([DS.dpath '\' DS.abfFn '.mat'],'file')
  [d,si]=matDload([DS.dpath '\' DS.abfFn '.mat'],'start',start,'stop',stop,'channels',rawChAnNm);
elseif exist([DS.dpath '\' DS.abfFn '.raw'],'file')
  rawFInfo=rawinfo('filename',[DS.dpath '\' DS.abfFn '.raw'],'print','no');
  d=rawload([DS.dpath '\' DS.abfFn '.raw'],rawChAnNm,[start stop]*1000,rawFInfo);
  % put into array and convert to mV
  d=cat(2,d{:})/1000;
  tmp=rawfi2genfi(rawFInfo,'channels',AP.rawChAnNm);
  si=tmp.si;
else
  [d,si]=abfload([DS.dpath '\' DS.abfFn '.abf'],'start',start,'stop',stop,'channels',rawChAnNm);
end
if DS.rawSignalInverted
  d=-1*d;
end
% downsmple and detrend (channel by channel for memory reasons) data 
% - many channels have offsets
d=d(1:sampFac:end,:);
si=si*sampFac;
for h=1:size(d,2)
  d(:,h)=d(:,h)-mean(d(:,h));
end

% start indices into d, also to be used as end indices (except for 
% last entry, which is an end index only)
intvPts=floor(linspace(1,size(d,1)+1,nRow+1));
% corresponding time points in real units (ms) and with the start of neural
% RECORDING (as opposed to local excerpt) as reference point
intv=WP.boe+discrete2cont(intvPts,si*.001,'intv',0);
intvDiff=intv(2)-intv(1);

for g=1:nRow
  subplot(nRow,1,g);
  % the combination of the following lines makes the subplots fill out 
  % pretty much the whole page
  tmpP=get(gca,'position');
  tmpPO=get(gca,'outerposition');
  % stretch x axis to the limits defined by outerposition...
  set(gca,'position',[tmpPO(1) tmpP(2) tmpPO(3) tmpP(4)])
  % but don't do this for y limits because the first and last subplots will
  % be distorted, and use rexy instead
  rexy('xfac',.95,'yfac',1.5);

  pllplot(d(intvPts(g):intvPts(g+1)-1,:),'si',si,'spacing','fixed','dy',dy,'ylim',yl,'noscb',1,'ylab','mV');

  % plot princ chan in dark blue
  c=flipud(get(gca,'children'));
  set(c(local_pcIx),'color',[0 0 .5]);
  
  % draw colored rectangles in background to indicate behavior
  % - x axis scaling by pllplot is arbitrary, depending on compression 
  % factor, so we have to map real time units of the time stamps to the x 
  % axis limits
  xl=get(gca,'xlim');
  % - multiply time stamps in ms by this factor to obtain the value for
  % current x axis
  xFac=diff(xl)/intvDiff;
  
  % - index to time stamps within current subplot limits 
  % + the last one of the previous subplot, if any. Not ethat the logic
  % below requires a time stamp at exactly t=0 in etsl
  tsIx=find(etsl(:,1)>=intv(g) & etsl(:,1)<intv(g+1));
  if ~isempty(tsIx)
    lastIx=tsIx(end);  
    if tsIx(1)>1
      tsIx=[tsIx(1)-1; tsIx];
    end
  else
    tsIx=lastIx;
  end

  for tsi=1:length(tsIx)
    col=briCol(etsl(tsIx(tsi),etslc.tagCol),:);
    % x limits of current patch in real units, replicated and arranged as 
    % 4-element column array to match ylArr
    xlArr=reshape(cumsum(etsl(tsIx([tsi tsi]),[etslc.tsCol etslc.durCol]),2),4,1);
    ph(tsi)=patch(xFac*(xlArr-intv(g)),ylArr,col);
    set(ph(tsi),'edgecolor','none');
  end
  
  % push rectangles to backgound
  c=get(gca,'children');
  set(gca,'children',circshift(c,nChan));

  if g==1,
    title([DS.abfFn ', ' DS.aName]);
  end
  if g==nRow
    set(gca,'FontSize',8);
    axis on
  end
end

if ~isempty(AP.printas{1}),
  for i=1:length(AP.printas)
    pa=AP.printas{i};
    if strfind(pa,'ps'), ext='.ps';
    elseif strfind(pa,'jpeg'),
      ext='.jpg';
      % override quality setting
      pa='-djpeg95';
    else
      ext='';
    end
  end
  print(pa,'-r400',[WP.figName '_' partJob ext]);
end

