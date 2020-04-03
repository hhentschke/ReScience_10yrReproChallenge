% This is a relatively quick-and-dirty routine which produces dose-response
% curves of the beta 3 knockin mice. For each given parameter, data are
% extracted for one electrode, both behaviors and both genotypes and are
% plotted versus etomidate concentration. Furthermore, linear regression
% plus F-tests are performed and the whole shabang gets plotted. The tests
% produce p-values (although confidence intervals would be better)
% informing you whether you may think that 
% - etomidate has an effect
% - wt and knockin differ under coontrol
% - wt and knockin differ in their sensitivity to etomidate
% and so on



% ---------- prelims ---------------------------------------------
labelscale('fontSz',16,'scaleFac',.8,'lineW',2.5,'markSz',16); 
curFigPath='c:\projects\hippocampus\rmouse\rawFig';
printas='-djpeg95';
printas='-dpsc2';
printas=[];


% use regstats to compute linear regression
% also use genotype as categorical var?
% glm to do the same?
% use non-normalized values & include control measurements??

% name of parameter of interest
% curRv='thgaeCCPeakMn_auto'
% curRv='thgaeCCPeakPhaseMn_auto'
% curRv='rawThPEMn_auto'
% curRv='rawGaPEMn_auto'
curRv='rawPPeakTMn_auto'


% recSite 0 is at fissure
recSite=0;
% recSite=-.40;
% loop over these files:
rfn={'rcomb_wtki_5','rcomb_wtki_10','rcomb_wtki_15'};
conc=[5 10 15]';


% -> here's an excerpt from combine_r.m explaning the structure of Rraw
% % struct holding collected results: R
% tmplt=cell(length(behav),length(rv),n3);
% tmplt(:)={[]};
% % all of the following fields are 3d cell arrays:
% % - row=behavior (in order listed above)
% % - col=parameter (in order listed above)
% % - slice=genotype (in order listed above)
% % - each element of the cell array contains this type of data:
% % collected data: 2d arr, 1st col electrode pos, 2nd col value (drug
% % exprmnts: 3rd+ cols = values ar var. concentrations)
% R.d=tmplt;       
wtIx=1;
kiIx=2;

% ---------- collect data ---------------------------------------------
% results variable:
mnD=[];
stD=[];
% individual control measurements 
ctrlMn=[];
ctrlSt=[];
% d is a cell array containing measurements from individual animals
% rows | columns: behavior | genotypye
% each cell contains an array (conc | dep var)
d=cell(2,2);
d(:)={[]};
for fi=1:length(rfn)
  load(rfn{fi});
  % loop over genotypes
  for gti=1:2
    % loop over behaviors
    for bi=1:2
      tmpd=Rraw.d{bi,strmatch(curRv,rv),gti};
      % pick rec Site
      tmpd(tmpd(:,1)~=recSite,:)=[];
      % control values without nans
      tmpd1=tmpd(isfinite(tmpd(:,2)),2);
      d{bi,gti}=cat(1,d{bi,gti},[repmat(0,size(tmpd1))  tmpd1]); 
      % drug values without nans
      tmpd1=tmpd(isfinite(tmpd(:,3)),3);
      d{bi,gti}=cat(1,d{bi,gti},[repmat(conc(fi),size(tmpd1))  tmpd1]);
      % compute means and std for plots
      ctrlMn(fi,gti,bi)=nanmean(tmpd(:,2));
      ctrlSt(fi,gti,bi)=nanstd(tmpd(:,2));
      tmpd=tmpd(:,3);
      mnD(fi,gti,bi)=nanmean(tmpd);
      stD(fi,gti,bi)=nanstd(tmpd);
    end
  end
end
clear tmp*

% ---------- fit data ---------------------------------------------
regressLineData=[];
for bi=1:2
  % --- first test (for each genotype separately): are slopes different from 0?
  % WT (adorn with ones)
  idepD=[ones(size(d{bi,wtIx}(:,1)))  d{bi,wtIx}(:,1)];
  % do it!
  [rcoeff,confI,tmpnix,tmpnix2,stats]=regress(d{bi,wtIx}(:,2),idepD);
  % put out a little information
  disp(['linear regression WT: r^2=' num2str(stats(1)) '; p=' num2str(stats(3))]);

  % MUT (adorn with ones)
  idepD=[ones(size(d{bi,kiIx}(:,1)))  d{bi,kiIx}(:,1)];
  % do it!
  [rcoeff,confI,tmpnix,tmpnix2,stats]=regress(d{bi,kiIx}(:,2),idepD);
  % put out a little information
  disp(['linear regression MU: r^2=' num2str(stats(1)) '; p=' num2str(stats(3))]);

  % ------ second & third test: slopes & offsets different between wt and ki?
  % --- PART I: unrestrained individual and combo fits 
  ft = fittype('a+x*b' ,...
    'dependent',{'y'},'independent',{'x'},...
    'coefficients',{'a','b'});
  % starting values for parameters
  st = [1 0];
  fo = fitoptions('method','NonlinearLeastSquares','Startpoint',st);
  % individual fits
  fit1=fit(d{bi,wtIx}(:,1),d{bi,wtIx}(:,2),ft,fo);
  fit2=fit(d{bi,kiIx}(:,1),d{bi,kiIx}(:,2),ft,fo);
  % combo fit
  fit12=fit(cat(1,d{bi,1}(:,1),d{bi,2}(:,1)),...
    cat(1,d{bi,1}(:,2),d{bi,2}(:,2)),ft,fo);


  % --- PART II: redo individual fits & F-test using fixed slope
  % of combo fit as value for null hypothesis
  tmpVal=coeffvalues(fit12);
  tmpVal=num2str(tmpVal(2));
  ft = fittype(['a+x*' tmpVal],...
    'dependent',{'y'},'independent',{'x'},...
    'coefficients',{'a'});
  % starting values for parameters
  st = [1];
  fo = fitoptions('method','NonlinearLeastSquares','Startpoint',st);
  fit1_fixSlope=fit(d{bi,wtIx}(:,1),d{bi,wtIx}(:,2),ft,fo);
  fit2_fixSlope=fit(d{bi,kiIx}(:,1),d{bi,kiIx}(:,2),ft,fo);
  
  % statistical test for similarity of slopes (cfit objects will be 'fevaluated' automatically)
  % ** note that we compare combinations of separately fitted data, hence
  % the degrees of freedom are # of data pts -2* # of free pars
  [p,F,radj1,radj2]=curvecomp2(...
    [[d{bi,wtIx}; d{bi,kiIx}]  [fit1_fixSlope(d{bi,wtIx}(:,1)); fit2_fixSlope(d{bi,kiIx}(:,1))]],...
    [[d{bi,wtIx}; d{bi,kiIx}]  [fit1(d{bi,wtIx}(:,1)); fit2(d{bi,kiIx}(:,1))]],...
    length(d{bi,wtIx}(:,1))+length(d{bi,kiIx}(:,1))-2*1,...
    length(d{bi,wtIx}(:,1))+length(d{bi,kiIx}(:,1))-2*2);
  
  urtext(['p=' num2str(p,'%1.3f')],.85,'fontsize',12);
  if p<.05
    disp(['***** H0 (identity of slopes) rejected, p= ' num2str(p)]);
    %       if p<.01, urtext('**',.9,'fontsize',20);
    %       else urtext('*',.9,'fontsize',20);
    %       end
  else
    disp(['H0 (identity of slopes) not rejected, p= ' num2str(p)]);
    %       urtext('n.s.');
  end

  % --- PART III: redo individual fits & F-test using fixed offset 
  % of combo fit as value for null hypothesis
  tmpVal=coeffvalues(fit12);
  tmpVal=num2str(tmpVal(1));
  ft = fittype([tmpVal '+x*b'],...
    'dependent',{'y'},'independent',{'x'},...
    'coefficients',{'b'});
  % starting values for parameters
  st = [0];
  fo = fitoptions('method','NonlinearLeastSquares','Startpoint',st);
  fit1_fixOffs=fit(d{bi,wtIx}(:,1),d{bi,wtIx}(:,2),ft,fo);
  fit2_fixOffs=fit(d{bi,kiIx}(:,1),d{bi,kiIx}(:,2),ft,fo);
  
  % statistical test for similarity of slopes (cfit objects will be 'fevaluated' automatically)
  % ** note that we compare combinations of separately fitted data, hence
  % the degrees of freedom are # of data pts -2* # of free pars
  [p,F,radj1,radj2]=curvecomp2(...
    [[d{bi,wtIx}; d{bi,kiIx}]  [fit1_fixOffs(d{bi,wtIx}(:,1)); fit2_fixOffs(d{bi,kiIx}(:,1))]],...
    [[d{bi,wtIx}; d{bi,kiIx}]  [fit1(d{bi,wtIx}(:,1)); fit2(d{bi,kiIx}(:,1))]],...
    length(d{bi,wtIx}(:,1))+length(d{bi,kiIx}(:,1))-2*1,...
    length(d{bi,wtIx}(:,1))+length(d{bi,kiIx}(:,1))-2*2);
  
  urtext(['p=' num2str(p,'%1.3f')],.85,'fontsize',12);
  if p<.05
    disp(['***** H0 (identity of offset) rejected, p= ' num2str(p)]);
    %       if p<.01, urtext('**',.9,'fontsize',20);
    %       else urtext('*',.9,'fontsize',20);
    %       end
  else
    disp(['H0 (identity of offset) not rejected, p= ' num2str(p)]);
    %       urtext('n.s.');
  end
  % collect some data for plots
  regressLineData(1:2,1:2,bi)=[fit1([0; conc(end)]) fit2([0; conc(end)])];

end


figure(1), clf, orient tall
% set(gcf,'pos',[659   382   651   719]);
set(gcf,'pos',[274    59   473   639]);
for bi=1:2
  subplot(2,1,bi), hold on
  ebh=errorbar(conc*[1 1],mnD(:,:,bi),stD(:,:,bi),'o');
  set(ebh(1),'color',[0 .5 .1],'markerfacecolor',[0 .5 .1]);
  set(ebh(2),'marker','^','color','b','markerfacecolor','b');
  if ~isempty(ctrlMn)
    ebh2=errorbar(conc*[0 0],ctrlMn(:,:,bi),ctrlSt(:,:,bi),'o');
    set(ebh2(1),'color',[0 .5 .1],'markerfacecolor',[0 .5 .1]);
    set(ebh2(2),'marker','^','color','b','markerfacecolor','b');
  end
  % now plot the regression line
  lh=line([0 conc(end)]',regressLineData(:,1,bi),'color',[0 .5 .1]);
  lh=line([0 conc(end)]',regressLineData(:,2,bi),'color','b');
  nicexyax;
  yl=get(gca,'ylim');
  % § accomodate negative ylims (thgaeCCPeakPhaseMn)
  if ~isempty(ctrlMn)
    set(gca,'xlim',[-1.5 16.5],'xtick',[0; conc],'ylim',[5*double(strcmpi(curRv,'rawPPeakTMn_auto')) yl(2)]);
  else
    set(gca,'xlim',[3.5 16.5],'xtick',conc,'ylim',[5*double(strcmpi(curRv,'rawPPeakTMn_auto')) yl(2)]);
  end
  rexy('ax',gca,'xfac',.75);
%   if bi==2, 
%     lh=legend(ebh,{'wt','{\beta}_3N265M'},'location','southwest');
%     % set(lh,'fontsize',14);
%   end
  xlabel('etomidate concentration (mg/kg)');
  ylabel([curRv ', ' RInfo(1).level{bi}]);

end
subpax(gcf)

if ~isempty(printas)
  print(printas,[curFigPath '\' curRv]);
end

% 'thgaeCCPeakMn_auto';...
% 'thgaeCCZScore_auto';...
% 'rawPPeakMn_auto';...
% 'rawPPeakTMn_auto';...
% 'rawThPEMn_auto';...
% 'rawThNarrowPEMn_auto';...
% 'rawBePEMn_auto';...
% 'rawGaPEMn_auto';...
% 'rawGaNarrowPEMn_auto';...
% 'rawRiPEMn_auto';...
% 'thPosPeakCvAMn_auto';...
% 'thNegPeakCvIPIMn_auto';...
% 'gaeThPEMn_auto';...
% 'gaeThNarrowPEMn_auto';...
% 'thCCPeakMn_cross';...
% 'thCCPeakTMn_cross';...
% 'rawCohMnTh_cross';...
% 'rawCohMnThNarrow_cross';...
% 'rawCohMnGa_cross';...
