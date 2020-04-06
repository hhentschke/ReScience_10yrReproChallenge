% An autonomous routine producing a handful of plots illustrating details of 
% crosscorrelation theta-gammaEnv & Matt's shuffle statistics. Needs output
% (=matfile) from seg_thgaeCC_demo_dgen

figdir='d:\projects\hippocampus\rmouse\paper_atropine\rawFig\';
% figdir='d:\hh\';
load([figdir 'seg_thgaeCC_demo_dgen']); 

% --- specific settings 
% number of shuffled versions of excerpt to display in 'raw' fig
nShuffP=3;
% excerpt (relative to segment borders) to plot in s
segIntv=[0 2]+1;
segIntvPts=cont2discrete(segIntv,si*1e-6);
segIntvPts=segIntvPts(1):segIntvPts(2);
% limit of freq for spectral plots
limFreq=[0 90];
% lags of cc to plot in ms
lagLim=500;
% 
ccHistBin=[.0:.02:.6];

% --- graphics
labelscale('fontSz',8,'scaleFac',.35,'lineW',1,'markSz',8); 
ornt='portrait';

printas='-dps2';
% printas=[];
figName=mfilename;

% 1. real data & shuffled data segments
% thin lines here
labelscale('fontSz',8,'scaleFac',.35,'lineW',.9,'markSz',8); 
% also, cut down traces
thetaD=thetaD(segIntvPts,:);
gammaEnvD=gammaEnvD(segIntvPts,:);
shGammaEnvD=shGammaEnvD(segIntvPts,:);
figure(1), clf,
pllplot([thetaD gammaEnvD shGammaEnvD(:,1:nShuffP)],'si',si);
if ~isempty(printas), print(printas,[figdir figName '_raw']); end

return

% 2. spectra
figure(2), clf, hold on
for g=1:nShuffP
  [P,F]=fspecp(shGammaEnvD(:,g),si,'limFreq',limFreq);
  ph(g)=plot(F,P);  
  set(ph(g),'color',[.4 .4 .4]);
end
[P,F]=fspecp(gammaEnvD,si,'limFreq',limFreq);
ph2=plot(F,P,'k');
set(ph,'linewidth',get(ph2,'linewidth')*.5);
axis tight
set(gca,'yscale','log','xscale','log')
set(gca,'xtick',[1 2.5 5 10 20 40 80]);
set(gca,'ytick',10.^[-10:0]);
if ~isempty(printas), print(printas,[figdir figName '_spec']); end

% 3. CC
figure(3), clf, hold on
lag=(size(sh_segCC,1)-1)/2;
lag=-lag:lag;
lag=lag*si*.001;
lagIx=lag>= -1*lagLim & lag<=lagLim;
ph=plot(lag(lagIx),sh_segCC(lagIx,2:nShuffP+1));
set(ph,'color',[.4 .4 .4]);
ph2=plot(lag(lagIx),sh_segCC(lagIx,1),'k');
set(ph,'linewidth',get(ph2,'linewidth')*.5);
set(gca,'ylim',max(abs(sh_segCC(lagIx,1)))*[-1.1 1.1]);
set(gca,'xtick',[-1500:250:1500]);
lh=line(get(gca,'xlim'),[0 0],'linestyle','--','linewidth',get(ph2,'linewidth')*.5,'color','k');
lh=line([0 0],get(gca,'ylim'),'linestyle','--','linewidth',get(ph2,'linewidth')*.5,'color','k');    
if ~isempty(printas), print(printas,[figdir figName '_CC']); end  

% 4. CC peaks
figure(4), clf, hold on
% determine peaks only on excerpts
tmpr=evdeal(sh_segCC(lagIx,:),'idx','minmaxpeak');
% spit out mean+std
disp(['mean shuffled peak CC: ' num2str(mean(tmpr.minPeak(2:end))*-1)]);
disp(['std shuffled peak CC: ' num2str(std(tmpr.minPeak(2:end)))]);
[n,x]=hist(tmpr.minPeak(2:end)*-1,ccHistBin);
bh=bar(x,n,1.0,'k');
niceyuax;
set(gca,'xlim',x([1 end]))
ph=plot(tmpr.minPeak(1)*-1,1,'kv');
if ~isempty(printas), print(printas,[figdir figName '_CCPeaks']); end  

% 5. cum hist all segments
figure(5), clf, hold on
ph=plot(shBins,shCumH);
set(ph,'color',[.4 .4 .4]);
ph=plot(realBins,realCumH,'k');
axis tight;
if ~isempty(printas), print(printas,[figdir figName '_CCPeakCumH']); end  

% 6. zscore
figure(6), clf, hold on
sigZ=2.5;
sigN=length(find(cccZScore>sigZ));
[n,x]=hist(cccZScore,[0:.25:ceil(max(cccZScore))+.25]);
bh=bar(x,n,1.0,'k');
set(gca,'xtick',[0:1:20]);
niceyuax;
plot(sigZ,10,'kv');
if ~isempty(printas), print(printas,[figdir figName '_ZScoreH']); end  


