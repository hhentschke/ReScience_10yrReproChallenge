% This script collects - into 5d-matrix D - selected results
% variables as computed by rmouse.m and combined for all sessions by 
% combine4b.m. D has a layout facilitating comparisons drug vs. control and
% correlations between the different results parameters.
% There are multiple options to plot the results. 
% To run this
% script, load the contents of auto.mat or cross.mat. Either of these mat 
% files contains a results variable R.
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

global D

% mm, limits of electrode depth (inclusive; slm=0, dorsal ones negative, ventral ones positive)
depthLim=[-.6 0];
nRecSites=7;

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
currRv={'thgaeCCPeakMn','rawThPEMn','rawGaPEMn'};
% currRv={'thgaeCCPeakMn','rawThPEMn'};
% currRv={'rawThNarrowPEMn','gaeThNarrowPEMn'};

[nada,rvIx]=intersect(rv,currRv);
% watch out - intersect sorts alphabetically
currRv=nada;
% rvIx=sort(rvIx); this is nonsense!
nCurrRv=length(currRv);

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

% --- collect data 
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


% --- data reduction and plot

plotType=2;

switch plotType
  case 3
    % electrodes: data points (symbol-coded)
    % rv: as many as specified above
    % behav: subplot
    % indv: in different plots
    % drug: all, color-coded
    % transformation required:
    % | drug | rv | el | behav | indv | -> | el | rv | indv |
    labelscale('fontSz',12,'scaleFac',1,'lineW',2,'markSz',8);
    figure(1), clf,
    for bInd=1:nBehav
      for drugInd=1:nDrug
        % following line below, we have as many columns as rv, 
        % all rec sites of individuals concatenated, individuals in slices
        tmpD=permute(D(drugInd,:,:,bInd,:),[3 2 5 1 4]);
        for iInd=1:nIndv
          indvD=tmpD(:,:,iInd);
          indvD(any(isnan(indvD),2),:)=[];
          [r,p]=corrcoef(indvD);
          currRv
          drugInd
          behav{bInd}
          p
          plot(indvD','o-')
          set(gca,'yscale','log');
          pause
        end
      end
    end
  case 2
    % electrodes: all, color-coded
    % rv: as many as specified above
    % behav: subplot
    % indv: data points
    % drug: all, symbol-coded
    % transformation required:
    % | drug | rv | el | behav | indv | -> | el | indv | rv |
    labelscale('fontSz',12,'scaleFac',1,'lineW',2,'markSz',8);
    figure(1), clf,
    for bInd=1:nBehav
      for drugInd=1:nDrug
        tmpD=permute(D(drugInd,:,:,bInd,:),[3 5 2 1 4]);
        % now we have as many columns as individuals, all rec sites of individuals concatenated, result pars in slices
        tmpD=reshape(tmpD,[nRecSites*nIndv,nCurrRv]);
        tmpD(any(isnan(tmpD),2),:)=[];
        % [r,p]=corrcoef(tmpD);
        currRv
        drugInd
        behav{bInd}
        % p
        plot(tmpD','o-')
        set(gca,'yscale','log');
        % **** ATTENTION: indices to tmpD arbitrary!
        [b,dev,stats]=glmfit(tmpD(:,[1 2]),tmpD(:,[3]),'normal');
        disp('\\\\\\\ generalized linear model: coefficients | standard errors | p-values \\\\\\\');
        kix=1;
        disp(['constant: ' sprintf(' %4.2f | %4.2f | %1.5f',[b(kix) stats.se(kix) stats.p(kix)])]);
        % one of the rv is the dependent var, therefore loop only up to nCurrRv
        for kix=2:nCurrRv
          disp([currRv{kix-1} ': ' sprintf(' %4.2f | %4.2f | %1.5f',[b(kix) stats.se(kix) stats.p(kix)])]);
        end
        pause
      end
    end
  case 1
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
end
