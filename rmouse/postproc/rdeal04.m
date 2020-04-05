% This script does a principal components analysis of selected results
% variables as computed by rmouse.m and combined for all sessions by 
% combine4b.m. To run this script, load the contents of auto.mat or cross.mat. 
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

global D

% mm, limits of electrode depth (inclusive; slm=0, dorsal ones negative, ventral ones positive)
depthLim=[-.6 0];
nRecSites=7;

depthLim=[ 0 0];
nRecSites=1;

% assume
% - all parameters were analyzed for each animal/conc/behav
uIndv=unique(R.indv{1,1});
nIndv=length(uIndv);
% number of drug conditions
nDrug=2;
% 
nBehav=length(behav);
% 
currRv={'thgaeCCPeakMn','rawThNarrowPEMn','gaeThNarrowPEMn'};
currRv={'rawDePEMn','rawThPEMn','rawBePEMn','rawGaPEMn'};
% currRv={'rawBePEMn','rawGaPEMn'};
[nada,rvIx]=intersect(rv,currRv);
% watch out - intersect sorts alphabetically
currRv=nada;
% rvIx=sort(rvIx); this is nonsense!
nCurrRv=length(currRv);
  
% number of observations=rows: N animals * 2 behaviors * 2 drug conditions
% number of characteristics=columns: number of selected electrodes * number of
% parameters
D=repmat(nan,nIndv*nDrug*nBehav,nRecSites*nCurrRv);

% colors and symbols for plotting
pset={'b','m';'o','s'};
% multiple small figs
labelscale('fontSz',8,'scaleFac',1,'lineW',2,'markSz',5);
for rvInd=1:nCurrRv
    for iInd=1:nIndv
      for drugInd=1:nDrug
        for bInd=1:nBehav
          rowIx=(iInd-1)*nDrug*nBehav + (drugInd-1)*nDrug + bInd;
          colIx=(rvInd-1)*nRecSites+1:rvInd*nRecSites;
          % tag for each of the four possible combos of behavior and drug cond
          tag(rowIx,1)=bInd+10*drugInd;

          % tag(rowIx,1)=drugInd;
          
          % identify data belonging to current session
          csInd=find(R.indv{bInd,rvIx(rvInd)}==uIndv(iInd));
          % 1st column rec depth, 2nd col ctrl, 3rd col drug
          % ** have rec depths in um because of that stupid rounding error
          recSites=floor(R.d{bInd,rvIx(rvInd)}(csInd,1)*1000);
          ds=R.d{bInd,rvIx(rvInd)}(csInd,drugInd+1);
          [nada,ix,ix2]=intersect(recSites,depthLim(1)*1000:100:depthLim(2)*1000);
          try
            % now put together
            D(rowIx,colIx(ix2))=ds(ix)';
          catch
            rvInd
            iInd
            drugInd
            bInd
          end

        end
      end
    end
end


PCexplore('tag',tag','omitMd','cha')