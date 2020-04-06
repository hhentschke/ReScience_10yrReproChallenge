% This script is a very much boiled-down version of rdeal.m,
% generating fits to & plots of thetaHi-gamma comodulation. 

% bar plots
labelscale('fontSz',6.5,'scaleFac',.34,'lineW',.7,'markSz',5); 
% this one for the somewhat larger lag plots
labelscale('fontSz',8,'scaleFac',.42,'lineW',1,'markSz',4.5); 

close all
fitArr={};
coll_ds12=[];

for nsi=1:2
  switch nsi
    case 1
      load immobile_thGaComod_auto.mat
    case 2
      load exploring_thGaComod_auto.mat
  end
  ds12=[ds1; ds2];

  % including data for x=0 makes sense for this parameter
  ds1ix=find(ds1(:,1)>-inf);
  ds2ix=find(ds2(:,1)>-inf);
  ds12ix=find(ds12(:,1)>-inf);
  % second-order polynomial
  ft_ = fittype('a + b*x + c*(x*x)' ,...
    'dependent',{'y'},'independent',{'x'},...
    'coefficients',{'a','b','c'});
  fo_ = fitoptions('method','NonlinearLeastSquares');
  % starting values for parameters
  st_ = [.15 -.3 1];

  % --- fit and put in fitArr
  set(fo_,'Startpoint',st_);
  [ds1f,gof1]=fit(ds1(ds1ix,1),ds1(ds1ix,2),ft_ ,fo_);
  [ds2f,gof2]=fit(ds2(ds2ix,1),ds2(ds2ix,2),ft_ ,fo_);
  [ds12f,gof12]=fit(ds12(ds12ix,1),ds12(ds12ix,2),ft_ ,fo_);

  % fill columns
  fitArr(nsi,1)={ds1f};
  fitArr(nsi,2)={ds2f};

  % --- part 4: create curves representing the fits
  % 1. all data pts & fit
  fitx=ds12(1,1):(ds12(end,1)-ds12(1,1))/200:ds12(end,1);
  % cfit objects will be 'fevaluated' automatically
  ds1fit=ds1f(fitx);
  ds2fit=ds2f(fitx);
  ds12fit=ds12f(fitx);

  % colors & symbols: there are as many groups of data as there are
  % levels of the second factor in plotColumnIx, so use these
  % levels' color specs
  pCol1='k';
  pSymb1='o';
  pCol2='b';
  pSymb2='s';


  [p,F,radj1,radj2]=curvecomp([ds1(ds1ix,:) ds1f(ds1(ds1ix,1))], ...
    [ds2(ds2ix,:) ds2f(ds2(ds2ix,1))],...
    [ds12(ds12ix,:) ds12f(ds12(ds12ix,1))],...
    length(ds1ix)-length(st_),length(ds2ix)-length(st_));

  if p<.05
    disp(['***** H0 (identity of depth-response profiles) rejected, p= ' num2str(p)]);
  else
    disp(['H0 (identity of depth-response profiles) not rejected, p= ' num2str(p)]);
    %       urtext('n.s.');
  end
  % goodness of fit
  [radj1 radj2]

  % paired ttests
  p=[];
  uRecSite=unique(ds12(:,1));
  nRecSite=length(uRecSite);
  for rsi=1:nRecSite
    ix1=ds1(:,1)==uRecSite(rsi);
    ix2=ds2(:,1)==uRecSite(rsi);
    if ~isequal(ix1,ix2)
      error('rec sites not same');
    end
    [h,p(rsi)]=ttest(ds1(ix1,2),ds2(ix2,2));
  end
  [uRecSite'; p*nRecSite]
  
  % rec site | drug | behav | subject
  coll_ds12=cat(1,coll_ds12,[ds1(:,1)  repmat(1,[size(ds1,1) 1])   repmat(nsi,[size(ds1,1) 1])   indv1   ds1(:,2)]);
  coll_ds12=cat(1,coll_ds12,[ds2(:,1)  repmat(2,[size(ds2,1) 1])   repmat(nsi,[size(ds2,1) 1])   indv2   ds2(:,2)]);

end
[rd,x]=regroup(coll_ds12);

compAxisPos=[.2 .2 .6 .6];
rdMn=nanmean(rd,4);
rdStd=nanstd(rd,0,4);

[axh,ph,lh]=errbarhh(x,rdMn(:,:,1),rdMn(:,:,2),rdStd(:,:,1),rdStd(:,:,2),...
  'pos',compAxisPos);
% colors of bars: there are as many groups of bars as there are
% levels of the second factor in plotColumnIx, so use these
% levels' color specs
set(ph(1,:),'facecolor',[.6 .6 .6],'barwidth',1.0);
set(ph(2,:),'facecolor','k','barwidth',1.0);
set(axh(:),'ytick',[],'XAxisLocation','top');

% fitArr=fliplr(fitArr);


for axhi=1:2
  subplot(axh(axhi)), hold on
  for curvei=1:2
    fitx=get(axh(axhi),'ylim');
    % cut fits a bit
    fitx=fitx+diff(fitx)*[.025 -.025];
    fitx=linspace(fitx(1),fitx(2),200);
    % Due to the way fitArr was filled we need to access it in
    % different ways depending on nMasterAnFacColIx
    curFit=fitArr{axhi,curvei};
    
    curFitY=curFit(fitx);
    fph=plot(curFitY,fitx,'k-');
    % again, if nMasterAnFacColIx==3 we want to do things a
    % little differently, namely have the two fits in each
    % double bar plot half window in different colors
    %       if 3
    %         set(fph,'color',...
    %           idv(plotColumnIx(2)).pCol{idv(plotColumnIx(2)).statsVal(curvei)});
    %       end
  end
  % if this is a bar plot shift BOTH fitted curves to backgound,
  % preserving their order
  set(gca,'chi',circshift(get(gca,'chi'),-2));
end

