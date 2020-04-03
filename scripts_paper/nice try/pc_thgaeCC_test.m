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



% --- II. 2D plot of thgaeCC (all animals, all electrodes) vs. new data (1st principal component)

% transformation required:
% | drug | rv | el | behav | indv | -> | el | indv | rv |

figure(11), clf, hold on
pstring={'co','bo','mv','rv'};
% pick behavior
for bInd=1:nBehav
  % pick drug condition 
  for drugInd=1:nDrug
    istring=[behav{bInd} ', drug condition: ' int2str(drugInd)];
    disp(istring);
    runInd=(bInd-1)*nBehav+drugInd;
    % tmpD shall have as many columns as individuals, all rec sites of
    % individuals concatenated, result pars in slices
    tmpD=permute(D(drugInd,:,:,bInd,:),[3 5 2 1 4]);
%     % now, normalize gamma and theta powers separately by respective powers at fissure
%     normFac=repmat(tmpD(princChInd,:,idepVarInd),[nRecSites 1 1]);
%     tmpD(:,:,idepVarInd)=tmpD(:,:,idepVarInd)./normFac;
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

    title(istring)
    labelscale('fontSz',8,'scaleFac',.6,'lineW',1.2,'markSz',6);
    
    ph=plot(pcD(:,1),pcD(:,2),pstring{runInd});
    
    
    
%     % princ components
%     [pcs,nd,vars,pcD]=PCexplore('data',pcD,'exploreMd',0,'normalize',0,'nPC',2);
%     % plot
%     figure(10+runInd), clf, hold on
%     title(istring)
%     labelscale('fontSz',8,'scaleFac',.3,'lineW',1.2,'markSz',6);
%     for ei=1:nRecSites
%       ix=find(recSTag==ei);
%       plot(nd(ix,1),tmpD(ix,3),['k' symb(ei)]);
%     end
    
    
  end
end
nicexyax

    
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
