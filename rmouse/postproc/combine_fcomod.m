function R=combine_fcomod
% extracts comodulograms from a collection of data sets (experiments) with 
% application of a drug, computes reduced data and and combines these. 
% needs ANPAR and DSET (=concatenation of individual and matching AP and DS)
% one experiment per column. Will mess up completely if row order of AP and 
% DS deviates from the following:
%       control
%       application
%       wash (optional)
%
% combine_fcomod works similar to the combineX (X any number from 1 to 4)
% functions in that it loops over data sets and collects data. The major
% difference is the nature of the data collected (single parameters in the
% other combine funcs, whole series here) and as a consequence, the
% structure of the final results variable R

global ANPAR DSET

% I. Set parameters ---------------------------------------------------------
% choose single results variable to collect and average/plot - 
% must be a field name of r (will be put in eval)
rv={'fComod'};
rvi=1;
% q MUST be auto
q='auto';
% choose behaviors to be treated (legal value of AP.segmentType)
behav={'immobile','exploring'};
behav={'immobile','exploring'};
nBehav=length(behav);

% frequency bands
AP_beta3_wtko;

% fBand=[AP.delta; [4 8]; [9 11]; AP.beta; [18 22]; AP.gamma];
% fbLabel={'\delta','\theta lo','\theta narrow','\beta','\beta narrow','\gamma'};

% fBand=[AP.delta; [4 8]; [8 11]; AP.beta; AP.gamma];
% fbLabel={'\delta','\theta lo','\theta narrow','\beta','\gamma'};

% fBand=[AP.theta; AP.gamma];
% fbLabel={'\theta','\gamma'};

fBand=[[4 8]; [8 11]; AP.gamma];
fbLabel={'\theta lo','\theta narrow','\gamma'};

% II. Various preparations ---------------------------------------------------------

rmouse_ini;

% --- plot stuff
% this generates variable segTypeGlobP containing global settings for plots of behaviors
stgp;
close all;
labelscale('fontSz',10,'scaleFac',1.0,'lineW',1.0,'markSz',6); 

% -------- PART III: collect information about recording sites & set up results var
[n1 n2 n3]=size(ANPAR);
[dn1 dn2 dn3]=size(DSET);
if ~isequal([n1 n2 n3],[dn1 dn2 dn3]), error('DSET and ANPAR must have equal dimension'); end
if n1>3, error(['current version of ' mfilename ' supports only a simple ''control | application (| recovery)''  scheme of drug experiments: ANPAR and DSET must have two or three lines']); end

% hp = 'helper parameters'
hp.elx=[];
% commChInd points to those LFP channels in dset which (i) were requested for analysis
% in original AP, and (ii) are available in ALL data sets of the current experiment
% (column of ANPAR and DSET)
% ** note again: the channels in all dsets of one experiment MUST be identical
commChInd=cell(1,n2);
for ci=1:n2
  commChNm=ANPAR(1,ci).rawChAnNm;
  for ri=1:n1
    AP=ANPAR(ri,ci);
    DS=DSET(ri,ci);
    rawCh=rmouse_chan;
    % hp.elx: union of all electrode positions in all experiments
    hp.elx=union(hp.elx,WP.elx);
    % let's suspect that for data won in one experiment
    % the number of LFP channels may be variable
    if ri==1
      % commChInd{ci}=1:nLFPCh;
      commChInd{ci}=AP.LFPccInd;
    else
      [commChInd{ci},ccni]=intersect(AP.LFPccInd,commChInd{ci});
      % [commChInd{ci},ccni]=intersect(1:nLFPCh,commChInd{ci});
      commChInd{ci}=AP.LFPccInd(sort(ccni));
    end
  end
end
hp.nCh=length(hp.elx);
% principal channel
hp.princChInd=find(hp.elx==0);
% final results field: linear cell array
f=cell(1,hp.nCh);

% ----- set up and/or list all major intermediate and final results variables ---
% 1. collected results: R
% - a struct array (as many elements as behaviors)
% - the field(s) are 1D cell arrays (one cell per channel); 
%   each data set will be embedded in the correct position
% - each cell is a 3D array:
%   - freq combo down columns, 
%   - column order control | drug | recovery (optional)
%   - different experiments in slices
f(:)={[]};
for i=1:nBehav
  R(i).r=f;
  R(i).indv=f;
end
 

% -------- PART IV: collection of data
% loop over data sets: 
% one experiment per column, concentration down the columns
% Some indices generated in this loop have to be saved for access of elements in the
% second (plotting) loop below. Assumption is that channels don't change within an
% experiment (are the same among behaviors and concentrations)
% Set up (preallocate, of sorts) the master index struct:
% 1. row & column indices into cell array
maIx(n2).rc=[];
% 2. slice index to 3D array within each cell (see further below)
maIx(n2).sl=[];
for ci=1:n2
  for ri=1:n1
    AP=ANPAR(ri,ci);
    DS=DSET(ri,ci);
    % need rmouse_chan for WP.elx and AP.LFPpcInd2
    rawCh=rmouse_chan;
    if ri==1
      % elix points to positions in final results cell array
      [elx,elix]=intersect(hp.elx,WP.elx(commChInd{ci}));
      % pre-generate position indices into master cell arr and current,
      % individual (slave) cell arr
      nLFPccCh=length(commChInd{ci});
      maIx(ci).rc=[ones(nLFPccCh,1) makecol(elix)];
      % ************************************************
      % in line below, note the difference to indexing in combine_rr  
      % (commented out below) due to the internal structure
      % of .fComod, which is different from all other spectral vars 
      % ************************************************      
      slave_sliceIx=(1:length(commChInd{ci}))';
      % slRC=[makecol(commChInd{ci}) makecol(commChInd{ci})];
      % slice index to 3D array within each cell
      maIx(ci).sl=repmat(nan,size(maIx(ci).rc,1),1);
    end
    % load results var..
    if isempty(strfind(AP.resPath,':')), AP.resPath=[WP.rootPath AP.resPath]; end    
    load([AP.resPath '\' AP.resFn],'r');
    % ..find behaviors..
    for bi=1:length(behav)
      % bix is index into r
      bix(bi)=strmatch(behav{bi},{r(:).segmentType});
    end
%     % define as narrow theta band the peak of theta at the principal channel
%     % during exploration +/- 1 Hz
%     currPeakF=r(strmatch('exploring',{r(:).segmentType})).rawPMnPeakT{AP.LFPpcInd2,AP.LFPpcInd2};
%     % ******* change if pos of theta narrow changes!
%     fBand(2,:)=currPeakF+[-1 1];    
    
    switch rv{1}
      case 'fComod'
        F=r(1).comF;
      otherwise
        error('Kamelle!')
    end
    % find freq indices
    disp('real freqs:')
    for fbi=1:size(fBand,1)
      fBandFIx{fbi}=find(r(1).comF>=fBand(fbi,1) & r(1).comF<=fBand(fbi,2));
      disp([num2str(F(fBandFIx{fbi}([1 end]))') ' Hz']);
    end
    [fComb,nFComb]=combin(size(fBand,1),'autoC',1);
    nanny=repmat(nan,nFComb,1);
    % and extract data
    for bi=1:length(bix)
      % tell what we're dealing with
      disp(['data: ' AP.resFn ', behavior: ' behav{bi} ', data:' rv{rvi}]);
      % tmpr is the original data
      eval(['tmpr=r(bix(bi)).' rv{rvi} ';']);
      allEmpty=isempty(tmpr);
      % loop over channels
      for rci=1:size(maIx(ci).rc,1)
        % slice index - identical for behaviors and concs
        if bi==1 & ri==1
          % funnily the size of [] in 3rd dim is 1, so catch that
          if isempty(R(bi).r{maIx(ci).rc(rci,1),maIx(ci).rc(rci,2)})
            maIx(ci).sl(rci)=1;
          else
            maIx(ci).sl(rci)=size(R(bi).r{maIx(ci).rc(rci,1),maIx(ci).rc(rci,2)},3)+1;
          end
        end
        % extract data, interpolate and put in proper position
        if allEmpty
          disp('**** no data, filling with nans');
          % this accounts for missing behaviors: fill with nans
          R(bi).r{maIx(ci).rc(rci,1),maIx(ci).rc(rci,2)}(:,ri,maIx(ci).sl(rci))=nanny;
        else
%           % theta deserves a special refinement: a strong narrowing. Within
%           % the standard theta freq band search for the peak of crosscorr
%           % values (theta x whole range). Do this for the principal channel only 
%           % and use this value for all others.
%           tmpMat=tmpr(fBandFIx{strmatch('\theta',fbLabel)},:,AP.LFPpcInd1);
%           [nada,peakF]=F(max(mean(tmpMat,2))+fBandFIx{strmatch('\theta',fbLabel)}(1)-1);
          fbArr=nanny;
          for fbi=1:nFComb
            % eliminate autocorrs becaue they are 1 and spoil 
            tmpr(:,:,slave_sliceIx(rci))=tmpr(:,:,slave_sliceIx(rci))...
              +diag(repmat(nan,size(tmpr,1),1),0);
            tmpMat=tmpr(fBandFIx{fComb(fbi,1)},fBandFIx{fComb(fbi,2)},slave_sliceIx(rci));
            fbArr(fbi)=mean(tmpMat(isfinite(tmpMat)));
%             if bi==2 & rci==12 & fbi==7,
%               figure(1), clf, imagesc(tmpr(:,:,slave_sliceIx(rci)),[-.3 .3]);
%               figure(2), clf, imagesc(tmpMat,[-.3 .3]);
%               disp('sdf');
%             end
          end
          R(bi).r{maIx(ci).rc(rci,1),maIx(ci).rc(rci,2)}(:,ri,maIx(ci).sl(rci))=fbArr;
          R(bi).indv{maIx(ci).rc(rci,1),maIx(ci).rc(rci,2)}(1,1,maIx(ci).sl(rci))=ci;
        end % if:allempty
      end % for:size(maIx(ci).rc,1)

    end % for:behav
  end % for:rows of ANPAR=concs
end % for:cols of ANPAR=experiments

save fComod_all R hp fBand fbLabel

