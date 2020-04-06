% produces time course contour plot of theta-gammaEnv CC of experimental
% session of choice
rmouse_ini;
AP_beta3_wtko;
figdir='d:\projects\hippocampus\rmouse\paper_atropine\rawFig\';

datdir=[WP.rootPath '\beta3_wtko\isoatr\wt2141_02523\thgaeCCexport\'];
figName=[mfilename '_isoAtr'];
cLim=[-.35 .35];

datdir=[WP.rootPath '\beta3_wtko\wt0002_04707\thgaeCCexport\'];
figName=[mfilename '_atr'];
cLim=[-.5 .5];


% --- specific settings 
% lags of cc to plot in ms
lagLim=250*[-1 1];

% --- graphics
labelscale('fontSz',10,'scaleFac',.8,'lineW',1,'markSz',8); 
ornt='portrait';

printas='-dps2';
% printas=[];

% list of all files
s=dir([datdir 'thgae*.mat']);
% global time
gT=[];
% global array of CC
gSegCC=[];
% loop over files, plot data on global time axis (which is inefficient, but
% collecting all data first and plotting them collectively is 
for k=1:length(s)
  load([datdir s(k).name]);
  % tmp_time is in terms of points; let's have the global time array in
  % points, too, plus rectime
  gT=[gT; tmp_time+cont2discrete(abfi.recTime(1),abfi.si/1e6,'intv',0)];
  gSegCC=[gSegCC segCC];
end

% cut down segCC according to lag range 

lag=-AP.ccLagPts:AP.ccLagPts;
lagIx=find(lag>=lagLim(1) & lag<=lagLim(2));
nLagIx=length(lagIx);
lag=lag(lagIx);
gSegCC=gSegCC(lagIx,:);


% determine contiguous blocks of data - since all segment-wise computed
% data do not fall in a fixed raster it is these blocks which have to be
% plotted individually
% - sort everything
[gT,tix]=sort(gT);
gSegCC=gSegCC(:,tix);
% - subtract time offset
gT=gT-gT(1);
% the elementary time step in the data in points
edt=(AP.ppSeg-AP.dftOlapPts);
% - time steps in gT
dt=diff(gT);
% - identify contiguous blocks
blockStart=[1; find(dt~=edt)+1];
blockEnd=[blockStart(2:end)-1;  length(gT)];
% - now convert time from pts to minutes
gT=discrete2cont(gT,abfi.si/6e7,'intv',0);
edt=discrete2cont(edt,abfi.si/6e7,'intv',0);
% - plot
figure(1), clf, orient(ornt); hold on

rexy('ax',gca,'xfac',.9,'yfac',.5);
colormap gray
for bi=1:length(blockStart)
  % account for the fact that the output of imagesc will be 1 unit wide 
  % in case a block  consists of one column of gSegCC only
  if blockStart(bi)==blockEnd(bi)
    x=gT(blockStart(bi))*[1 1]+[0 edt];
    z=[gSegCC(:,blockStart(bi)) repmat(cLim(2),nLagIx,1)];
  else
    x=gT([blockStart(bi) blockEnd(bi)]);
    z=gSegCC(:,blockStart(bi):blockEnd(bi));
  end
  imagesc(x,lag,z,cLim);
end


axis tight
set(gca,'ytick',-400:100:400,'xtick',0:20:300);
xlabel('Time (min)');
ylabel('Lag (ms)');
cbh=colorbar;
set(cbh,'ytick',-0.5:0.25:0.5);

if ~isempty(printas), print(printas,[figdir figName]); end

