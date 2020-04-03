% produces polar plots of CC (whole waveforms) generated
% from short segments (the ..iState files)
% - needs 

global WP

% which session?
expInd=1;
% downsampling factor
sampFac=4;

% collect DS and AP
collect_wtAtropine_av; 

% load data
D=[];
for condInd=1:2
  AP=ANPAR(condInd,expInd);
  DS=DSET(condInd,expInd);  
  load([WP.rootPath AP.resPath '\thgae_' DS.abfFn '_' AP.rawChPrincNm{:} '_exploring']);
  % load([WP.rootPath AP.resPath '\thgae_' DS.abfFn '_IN 8_exploring']);
  nSeg(condInd)=size(segCC,2);
  D=[D segCC];
end
AP_beta3_wtko_iStateTheta;

% D=D(41:end-40,:);
% AP.ccLagPts=AP.ccLagPts-41;

% make tag
tag=[zeros(nSeg(1),1);ones(nSeg(2),1)];
[n1,n2]=size(D);


osi=1000;
% time (ms)
t=discrete2cont([1:n1]'-AP.ccLagPts,osi*.001,'intv',0);
% time in radians for fitting purposes
tRad=t*.001*2*pi;

% determine peaks so we have good starting pars
r=evdeal(D,'idx','minmaxpeak');
rhStart=r.minPeak';
phStart=discrete2cont(r.minPeakT'-AP.ccLagPts,osi*.001,'intv',0)*.001*2*pi;

% cosine*gaussian
% ft=fittype('(rh*cos(om*tRad+ph)) * exp(-b*(tRad.^2))' ,...
%   'dependent',{'y'},'independent',{'tRad'},...
%   'coefficients',{'rh','om','ph','b'});

ft=fittype('(rh*cos(om*tRad+ph)) * exp(-b*((tRad+a).^2))' ,...
  'dependent',{'y'},'independent',{'tRad'},...
  'coefficients',{'rh','om','ph','b','a'});

fo=fitoptions('method','NonlinearLeastSquares');
% fo=fitoptions('method','NonlinearLeastSquares',...
% 'Lower',[-1 4 -10 eps],'Upper',[1 20 10 10]);

% array of starting values for parameters
% startPar = [rhStart repmat(9,n2,1) phStart repmat(.4,n2,1)];
startPar = [rhStart repmat(9,n2,1) phStart repmat(.4,n2,1) repmat(0,n2,1)];

coeff=repmat(nan,n2,size(startPar,2));
for segInd=1:500% n2
  y=D(:,segInd);
  set(fo,'Startpoint',startPar(segInd,:));
  yfit=fit(tRad,y,ft,fo);
  fifi=yfit(tRad);
  coeff(segInd,:)=coeffvalues(yfit);
%   clf,
%   plot(tRad,D(:,segInd));
%   hold on,
%   plot(tRad,fifi,'r');
%   niceyax;
%   pause;
end

% angle=lag: ph as computed above normalized by freq of sine
polar(coeff(:,3).*coeff(:,2),coeff(:,1),'o');

return



figure(44), clf
orient landscape
labelscale('fontSz',8,'scaleFac',1,'lineW',.5,'markSz',2); 
subplot(3,3,5), hold on, box on
% 2nd vs. 1st, control 
ph=plot(nd(tag==0,1),nd(tag==0,2),'ko');
% 2nd vs. 1st, atropine
ph=plot(nd(tag==1,1),nd(tag==1,2),'ko');
set(ph,'markerfacecolor','k');
nicexyax;
xl=get(gca,'xlim');
yl=get(gca,'ylim');
% --- 1d-histograms
bin1=linspace(xl(1),xl(2),nBin);
bin2=linspace(yl(1),yl(2),nBin);
% - 1st component (abscissa)
subplot(3,3,2), hold on
% control
[n,h]=hist(nd(tag==0,1),bin1);
% ** normalize for number of events
n=n/nSeg(1);
plot(h,n,'g-');
% atropine
[n,h]=hist(nd(tag==1,1),bin1);
% ** normalize for number of events
n=n/nSeg(2);
plot(h,n,'r-');
niceyuax;
set(gca,'xlim',xl);
set(gca,'xtick',[]);

% - 2nd component (ordinate)
subplot(3,3,6), hold on
% control
[n,h]=hist(nd(tag==0,2),bin2);
% ** normalize for number of events
n=n/nSeg(1);
plot(n,h,'g-');
% atropine
[n,h]=hist(nd(tag==1,2),bin2);
% ** normalize for number of events
n=n/nSeg(2);
plot(n,h,'r-');
nicexax;
tmpxl=get(gca,'xlim');
set(gca,'xlim',[0 tmpxl(2)]);
set(gca,'ytick',[]);


% --- difference of (normalized) densities
subplot(3,3,3)
colormap(coma('redblue'));
% ctrl
h_ctrl=hist2d(nd(tag==0,[1 2]),bin1,bin2)/nSeg(1);
h_atr=hist2d(nd(tag==1,[1 2]),bin1,bin2)/nSeg(2);
h=h_ctrl-h_atr;
ph=image(bin1,bin2,flipud(rot90(h)));
set(gca,'ydir','normal');
set(ph,'cdatamapping','scaled');
colorbar
% § set axis to values of scatter plot
