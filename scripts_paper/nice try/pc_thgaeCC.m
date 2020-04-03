% This script collects - into 5d-matrix D - selected results
% variables as computed by rmouse.m and combined for all sessions by 
% combine4b.m. D has a layout facilitating comparisons drug vs. control and
% correlations between the different results parameters.
% There are multiple options to plot the results. 
% Currently, the scripts generates 2 plots:
% i) 3D stem plot of CC theta-gammaEnv vs. theta power and gamma power (one
% individual)
% ii) 2D scatter plot of CC theta-gammaEnv vs. 1st principal component of
% [gamma power, theta power] (all animals)

% To run this script, load the contents of auto.mat or cross.mat. 
% Either of these mat files contains a results variable R.
% R is the collection of results from several experiments, won by running
% combine4b.m. It is a struct with (currently) 6 fields. Each field is a n-D cell
% array:
% - row=behavior (order as in behav)
% - col=parameter (order as in rv)
% - slice=genotype (1st Wt, 2nd mutant)
% Description of the fields (i.e., each individual cell's content):
%         d: collected data: 2d arr, 1st col electrode pos, 2nd col value (drug
%            exprmnts: 3rd+ cols = values ar var. concentrations)
%      indv: code for individual animal/session
%      ueix: 1d cell array; for each of the ue, these are the indices into the
%            corresponding R.d
%        ga: 'grand average': 2d arr, holding  mean|std|N
%     bstat:
%        ue:


printas='-dpsc2';
printas=[];
figdir='d:\projects\rmouse\paper_atropine\rawFig\';


% mm, limits of electrode depth (inclusive; slm=0, dorsal ones negative, ventral ones positive)
depthLim=[-.6 0];
princChInd=7;
nRecSites=7;
% symbols for rec sites (alveus first element, HF last)
symb=['o','o','^','s','s','*','*'];


% assume
% - all parameters were analyzed for each animal/conc/behav
uIndv=unique(R.indv{1,1});
nIndv=length(uIndv);
% number of drug conditions
nDrug=2;
% 
nBehav=length(behav);
% 
% currRv={'rawGaPEMn','gaeThNarrowPEMn'};
% currRv={'rawThPEMn','rawThNarrowPEMn'};
% currRv={'rawGaPEMn','gaeThPEMn'};
% currRv={'rawGaPEMn','rawThNarrowPEMn'};
% currRv={'rawGaPEMn','rawThPEMn'};
% currRv={'thCCPosPeakDecayMn','gaeThNarrowPEMn'};
% currRv={'thCCPosPeakDecayMn','rawThNarrowPEMn'};
currRv={'thgaeCCPeakMn','gaeThNarrowPEMn'};
currRv={'thgaeCCPeakMn','rawThNarrowPEMn'};
currRv={'thgaeCCPeakMn','rawThNarrowPEMn','gaeThNarrowPEMn'};
currRv={'thgaeCCPeakMn','rawThNarrowPEMn','rawGaPEMn'};
currRv={'thgaeCCZScore','rawThPEMn','rawGaPEMn'};
% currRv={'thgaeCCPeakMn','rawThPEMn','rawGaPEMn'};
% currRv={'thgaeCCPeakMn','rawThPEMn'};
% currRv={'rawThNarrowPEMn','gaeThNarrowPEMn'};

[nada,rvIx]=intersect(rv,currRv);
% watch out - intersect sorts alphabetically
currRv=nada;
nCurrRv=length(currRv);

% dependent var for regression and its index in currRv
depVar='thgaeCCZScore';
% depVar='thgaeCCPeakMn';
depVarInd=strmatch(depVar,currRv);
idepVarInd=setdiff(1:nCurrRv,depVarInd);

% ***********************************************************
% 1st row control
% 2nd row drug
% 1st col (abscissa): first par
% 2nd col (ordinate): second par
% along slices: electrodes
% along 4th dim: behavior
% along 5th dim: animals
D=repmat(nan,[nDrug nCurrRv nRecSites nBehav nIndv]);
% ***********************************************************

% --- 0. collect data 
for rvInd=1:nCurrRv
  for iInd=1:nIndv
    for drugInd=1:nDrug
      for bInd=1:nBehav
        % identify data belonging to current session
        csInd=find(R.indv{bInd,rvIx(rvInd)}==uIndv(iInd));
        % 1st column rec depth, 2nd col ctrl, 3rd col drug
        % ** have rec depths in um because of that stupid rounding error
        recSites=floor(R.d{bInd,rvIx(rvInd)}(csInd,1)*1000);
        ds=R.d{bInd,rvIx(rvInd)}(csInd,drugInd+1);
        [nada,ix,ix2]=intersect(recSites,depthLim(1)*1000:100:depthLim(2)*1000);
        % now put together
        D(drugInd,rvInd,ix2,bInd,iInd)=permute(ds(ix),[3 2 1]);
      end
    end
  end
end


% --- I. 3D stem plot of one individual, regression theta pow vs gamma pow
labelscale('fontSz',8,'scaleFac',.3,'lineW',1.8,'markSz',10);
figure(1), clf, 
% pick behavior
for bInd=1% 1:nBehav
  % pick drug condition 
  for drugInd=2 %:nDrug
    istring=[behav{bInd} ', drug condition: ' int2str(drugInd)];
    disp(istring);
    % following line below, we have as many columns as rv,
    % all rec sites of individuals concatenated, individuals in slices
    tmpD=permute(D(drugInd,:,:,bInd,:),[3 2 5 1 4]);
    % pick individual animal
  % for iInd=3
    for iInd=1:nIndv
      indvD=tmpD(:,:,iInd);
      indvD(any(isnan(indvD),2),:)=[];

      % ---- regression gamma vs theta
      % for this regression, gamma power is indep var and theta power dep var
      idepVar=[ones(size(indvD,1),1) indvD(:,idepVarInd(1))];
      [rcoeff,confI,nada1,nada2,stats]=regress(indvD(:,idepVarInd(2)),idepVar);
      % put out a little information
      disp(['linear regression ' currRv{idepVarInd(2)} '=B*' currRv{idepVarInd(1)}  '+C: r^2=' num2str(stats(1)) '; p=' num2str(stats(3))]);

      % regression thgaeCCZscore vs. PC1
      % var to be put into PC analysis shall contain only theta and gamma power
      pcD=indvD(:,idepVarInd);
      % princ components
      [pcs,nd,vars,pcD]=PCexplore('data',pcD,'exploreMd',0,'normalize',1,'nPC',2);
      % regression, thCCZScore is dep var and PC1 indep var
      idepVar=[ones(size(indvD,1),1) nd(:,1)];
      [rcoeff,confI,nada1,nada2,stats]=regress(indvD(:,depVarInd),idepVar);
      % put out a little information
      disp(['linear regression ' currRv{depVarInd} '=B*PC1+C: r^2=' num2str(stats(1)) '; p=' num2str(stats(3))]);
      pause

      % stem plot without markers
      stem3(indvD(:,idepVarInd(1)),indvD(:,idepVarInd(2)),indvD(:,depVarInd),'.k-');
      % zlabel(currRv{depVarInd});
      nicexyzax(15);
      % plot regression line in xy plane
      xl=get(gca,'xlim');
      yl=get(gca,'ylim');
      zl=get(gca,'zlim');      
      set(gca,'zlim',[0 zl(2)]);
      lh=line(xl,rcoeff(1)+rcoeff(2)*xl,'color',[.5 .5 .5],'linestyle','--');
      hold on
      % now plot symbols 
      for ei=1:size(indvD,1)
        stem3(indvD(ei,idepVarInd(1)),indvD(ei,idepVarInd(2)),indvD(ei,depVarInd),['k-' symb(ei)]);
      end
      
      set(gca,'CameraPosition',[-0.432767 -0.761979 4.03129],...
        'CameraTarget',[0.0341855 0.125353 0.518984]);

     
      % *** suspend glm because the distribution of theta-gammaEnvCC Z scores is very nonstandard (kind of U-shaped)
%       % glm
%       [b,dev,stats]=glmfit(indvD(:,setdiff(1:nCurrRv,depVarInd)),indvD(:,depVarInd),'normal');
%       disp('\\\\\\\ generalized linear model: coefficients | standard errors | p-values \\\\\\\');
%       kix=1;
%       disp(['constant: ' sprintf(' %4.2f | %4.2f | %1.5f',[b(kix) stats.se(kix) stats.p(kix)])]);
%       % one of the rv is the dependent var, therefore loop only up to nCurrRv
%       for kix=2:nCurrRv
%         disp([currRv{kix-1} ': ' sprintf(' %4.2f | %4.2f | %1.5f',[b(kix) stats.se(kix) stats.p(kix)])]);
%       end
      disp('')
    end
  end
end

if ~isempty(printas),
  print(printas,[figdir '3D_thetaGammaCorr']);
end

return
% --- II. 2D plot of thgaeCC (all animals, all electrodes) vs. new data (1st principal component)

% transformation required:
% | drug | rv | el | behav | indv | -> | el | indv | rv |
% pick behavior
for bInd=1:nBehav
  % pick drug condition 
  for drugInd=1:nDrug
    istring=[behav{bInd} ', drug condition: ' int2str(drugInd)];
    disp(istring);
    % tmpD shall have as many columns as individuals, all rec sites of
    % individuals concatenated, result pars in slices
    tmpD=permute(D(drugInd,:,:,bInd,:),[3 5 2 1 4]);
    % now, normalize gamma and theta powers separately by respective powers at fissure
    normFac=repmat(tmpD(princChInd,:,idepVarInd),[nRecSites 1 1]);
    tmpD(:,:,idepVarInd)=tmpD(:,:,idepVarInd)./normFac;
    % reshape such that individuals are concatenated
    tmpD=reshape(tmpD,[nRecSites*nIndv,nCurrRv]);
    % numerical tag for rec sites
    recSTag=repmat([1:nRecSites]',nIndv,1);
    % kill nans
    nanix=any(isnan(tmpD),2);
    tmpD(nanix,:)=[];
    recSTag(nanix)=[];
    % var to be put into PC analysis shall contain only theta and gamma power
    pcD=tmpD(:,idepVarInd);
    
    % princ components
    [pcs,nd,vars,pcD]=PCexplore('data',pcD,'exploreMd',0,'normalize',0,'nPC',2);
    
    % for this regression, PC1 is indep var and thgaeCCZScore dep var
    idepVar=[ones(size(nd,1),1) nd(:,1)];
    [rcoeff,confI,nada1,nada2,stats]=regress(tmpD(:,depVarInd),idepVar);
    % put out a little information
    disp(['linear regression ' currRv{depVarInd} '=B*<th&ga pow>+C: r^2=' num2str(stats(1)) '; p=' num2str(stats(3))]);

    
    
    % plot
    figure(10+(bInd-1)*nBehav+drugInd), clf, hold on
    title(istring)
    labelscale('fontSz',8,'scaleFac',.3,'lineW',1.2,'markSz',6);
    for ei=1:nRecSites
      ix=find(recSTag==ei);
      plot(nd(ix,1),tmpD(ix,3),['k' symb(ei)]);
    end
    nicexyax
    
    
    %     stem3(tmpD(:,idepVarInd(1)),tmpD(:,idepVarInd(2)),tmpD(:,depVarInd),'filled');
    %     [b,dev,stats]=glmfit(tmpD(:,setdiff(1:nCurrRv,depVarInd)),tmpD(:,depVarInd),'normal');
    %     disp('\\\\\\\ generalized linear model: coefficients | standard errors | p-values \\\\\\\');
    %     kix=1;
    %     disp(['constant: ' sprintf(' %4.2f | %4.2f | %1.5f',[b(kix) stats.se(kix) stats.p(kix)])]);
    %     % one of the rv is the dependent var, therefore loop only up to nCurrRv
    %     for kix=2:nCurrRv
    %       disp([currRv{kix-1} ': ' sprintf(' %4.2f | %4.2f | %1.5f',[b(kix) stats.se(kix) stats.p(kix)])]);
    %     end

  end
end

if ~isempty(printas),
  print(printas,[figdir 'thetaGammaCorr_PC']);
end


return
% electrode: princ only
% rv: first two as axes
% behav: expl & imm in separate subplots
% indv: data points
% drug: data points, color-coded
labelscale('fontSz',12,'scaleFac',1,'lineW',2,'markSz',8);
figure(1), clf,
for bInd=1:nBehav
  subplot(nBehav,nBehav,bInd), hold on
  for iInd=1:nIndv
    % princ chan only
    for recInd=7
      if 0
        % nrmlzd
        plot(D(:,1,recInd,bInd,iInd)./repmat(D(1,1,recInd,bInd,iInd),nDrug,1),...
          D(:,2,recInd,bInd,iInd)./repmat(D(1,2,recInd,bInd,iInd),nDrug,1),'o-');
      else
        ph=plot(D(:,1,recInd,bInd,iInd),...
          D(:,2,recInd,bInd,iInd),'o-');
        ph=plot(D(nDrug,1,recInd,bInd,iInd),...
          D(nDrug,2,recInd,bInd,iInd),'o-');
        set(ph,'markerfacecolor','b');
      end
    end
  end
  nicexyax;
  llim=[min(get(gca,'xlim'),get(gca,'ylim'))  max(get(gca,'xlim'),get(gca,'ylim'))];
  lh=line(llim,llim,'linestyle','--','color','r');
  xlabel(currRv{1});
  ylabel(currRv{2});
end
