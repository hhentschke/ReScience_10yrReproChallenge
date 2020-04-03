function R=combine_rr
% extracts and combines 'raw results' (rr) like spectra or crosscorrelations 
% from a collection of data sets (experiments) with application of a drug.
% needs ANPAR and DSET (=concatenation of individual and matching AP and DS)
% one experiment per column. Will mess up completely if row order of AP and 
% DS deviates from the following:
%       control
%       application
%       wash (optional)
%
% combine_rr works similar to the combineX (X any number from 1 to 4) functions
% in that it loops over data sets and collects data. The major difference is the nature of
% the data collected (single parameters in the other combine funcs, whole series here) and
% as a consequence, the structure of the final results variable R

global ANPAR DSET AP DS WP

% I. Set parameters ---------------------------------------------------------
% choose single results variable to collect and average/plot - 
% must be a field name of r (will be put in eval)
% currently, works only with power spectra (rawPMn and gaePMn)
rv={'gaePMn'};
rv={'rawPMn'};
% choose auto- or cross-channel results (cross not computed for theta, gamma env corr)
q='cross';
q='auto';
% choose behaviors to be treated (legal value of AP.segmentType)
behav={'immobile','exploring'};
nBehav=length(behav);

% spectra etc. will be up/downsampled to this resolution (Hz)
% -> needed because sampling intervals were not identical in all cases
fResol=.25;
% desired frequency range of spectra for plots
rfreq=[4:fResol:12];
rfreq=[1:fResol:95];
rfreq=[40:fResol:100];
rfreq=[100:fResol:250];
rfreq=[1:fResol:20];
rfreq=[20:fResol:100];
% frequency range within which to compute integral, peak, etc
% !must be contained in rfreq!
anfreq=rfreq;
% anfreq=[60:fResol:70];
% frequency range which shall be considered for normalization 
% !must be contained in rfreq!
normrfreq=rfreq; 
% kind of normalization of (difference) spectra:
% a) 'none'
% b) 'ctrlPow': for each channel, divide spectra by total spectral power in control case 
%    (= divide by integral (=area) of control spectrum)
%    -> this normalization emphasizes RELATIVE differences between control and drug
%    condition, even if the absolute differences are small
% c) 'diffPow': for each channel, normalize DIFFERENCE of spectra such that 
%    the difference (drug-ctrl) is equally 1
%    -> this normalization upscales channels with low absolute power

normType='diffPow';
normType='ctrlPow';
normType='none';
% shall log of power density be plotted?
logp=0;
% print?
printas=[];'-djpeg90';

% II. Various preparations ---------------------------------------------------------
rfreqix=1:length(rfreq);
% ** note that the two indices below are relative (local) indices
normrfreqix=find(rfreq>=normrfreq(1) & rfreq<=normrfreq(end));
anfreqix=find(rfreq>=anfreq(1) & rfreq<=anfreq(end));

rmouse_ini;

% --- plot stuff
% this generates variable segTypeGlobP containing global settings for plots of behaviors
stgp;
% in plots, set vertical lines at these frequency values
vlines=[4 8 12];
close all;
labelscale('fontSz',10,'scaleFac',1.0,'lineW',1.0,'markSz',6); 
% current plot scheme: two subplots per behavior, all in one single row,
% the first plot of each set control and drug together, the second difference
nPCols=nBehav*2;
dy=cell(1,nBehav);

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
      commChInd{ci}=AP.LFPccInd;
    else
      [commChInd{ci},ccni]=intersect(AP.LFPccInd,commChInd{ci});
      commChInd{ci}=AP.LFPccInd(sort(ccni));
    end
  end
end
hp.nCh=length(hp.elx);
% principal channel
hp.princChInd=find(hp.elx==0);
% final results field: linear or 2D cell array
switch q
  case 'auto'
    f=cell(1,hp.nCh);
  case 'cross'
    f=cell(hp.nCh);    
end

% ----- set up and/or list all major intermediate and final results variables ---
% 1. collected results: R
% - a struct array (as many elements as behaviors)
% - the field(s) are either
%   (i) 1D cell arrays if only auto-results are requested
%   (ii) 2D cell arrays if cross-results are requested (all non-redundant 
%   combinations of channels). Rows and columns span all recording sites
%   of all data sets; each data set will be embedded in the correct position
% - each cell is a 3D array:
%   - freq or time down columns, 
%   - column order control | drug | recovery (optional)
%   - different experiments in slices
f(:)={[]};
for i=1:nBehav
  R(i).r=f;
end
% 2. same, but each cell contains statistics of difference spectra obtained from
% all experiments: 
% - freq down columns
% - columns order:  mean | std 
% - slice order: drug-control | recov-control
mnRdiff=R;

% temporary var for plotting and derived results calculation
% indices: freq | channel | conc
plotR=[];

% 4. another helper variable containing difference spectra for each individual
% experiment: 
%   - freq down columns, 
%   - each channel (combo) one column 
%   - slice order: drug-ctrl | wash-ctrl
Rdiff=[];


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
    if ri==1
      % need rmouse_chan only for WP.elx
      rawCh=rmouse_chan;
      % elix points to positions in final results cell array
      [elx,elix]=intersect(hp.elx,WP.elx(commChInd{ci}));
      % pre-generate position indices into master cell arr and current,
      % individual (slave) cell arr
      nLFPccCh=length(commChInd{ci});
      switch q
        case 'auto'
          maIx(ci).rc=[ones(nLFPccCh,1) makecol(elix)];
          slRC=[makecol(commChInd{ci}) makecol(commChInd{ci})];
        case 'cross'
          % row & column index in 2D CC LFP cell arrays accessing upper
          % triangular part including main diagonal
          co=combin(nLFPccCh,'autoC',1);
          maIx(ci).rc=[elix(co(:,1)) elix(co(:,2))];
          slRC=[(commChInd{ci}(co(:,1)))' (commChInd{ci}(co(:,2)))'];
      end
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
    switch rv{1}
      case 'rawPMn'
        F=r(1).F;
      case 'gaePMn';
        F=r(1).gaeF;        
    end
    % find freqs of interest (end values must be a little wider because interpolation
    % will yield nans for x values outside interval)
    frix=find(F>=rfreq(1)-fResol & F<=rfreq(end)+fResol);
    % preparations for cutting out line hum (from complete spectrum)
    humF=60:60:F(end);
    for h=1:length(humF)
      % indices to immediately adjacent bins
      [nada,tmpix]=min(abs(F-humF(h)));
      % the ones to be replaced
      ix1a=tmpix-1:tmpix+1;
      ix1a(ix1a<1 | ix1a>length(F))=[];
      ix1{h}=ix1a;
      % the ones to compute mean from & replace with 
      ix2a=[tmpix-5:tmpix-2 tmpix+2:tmpix+5];
      ix2a(ix2a<1 | ix2a>length(F))=[];            
      ix2{h}=ix2a;
    end
    % and extract data
    for bi=1:length(bix)
      for rvi=1:length(rv)
        % tell what we're dealing with 
        disp(['data: ' AP.resFn ', behavior: ' behav{bi} ', data:' rv{rvi}]);      
        % tmpr is the original data 
        eval(['tmpr=r(bix(bi)).' rv{rvi} ';']);
        allEmpty=isequal(size(tmpr),[1 1]);
        % missing behaviors
        if allEmpty
          disp('**** no data, filling with nans');
          nanny=repmat(nan,rfreqix(end),1);
        end
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
          % this is just a final check
          if ~allEmpty
            if isempty(tmpr{slRC(rci,1),slRC(rci,2)}) | isnan(tmpr{slRC(rci,1),slRC(rci,2)})
              error('indexing error or faulty data');
            end
          end
          % extract data, interpolate and put in proper position
          if allEmpty 
            % this accounts for missing behaviors: fill with nans
            R(bi).r{maIx(ci).rc(rci,1),maIx(ci).rc(rci,2)}(rfreqix,ri,maIx(ci).sl(rci))=nanny;
          else
            % whole spectrum of current channel pair
            % ** note: for cross measures, compute absolute value here 
            tmpP=abs(tmpr{slRC(rci,1),slRC(rci,2)});
            for h=1:length(ix1)
              % this substitutes 60 Hz peaks and harmonics by mean of neighboring bins
              tmpP(ix1{h},:)=repmat(mean(tmpP(ix2{h},:),1),[3 1]);          
            end
            R(bi).r{maIx(ci).rc(rci,1),maIx(ci).rc(rci,2)}(rfreqix,ri,maIx(ci).sl(rci))=...
              interp1(F(frix),tmpP(frix),rfreq);
          end
        end % for:size(maIx(ci).rc,1)
      end % for:par
    end % for:behav
  end % for:rows of ANPAR=concs
end % for:cols of ANPAR=experiments


% -------- PART V: plot & precompute and collect normalized spectra and differences for average
% 1st row immobile, 2nd row exploring
% experimental sessions along columns
% 1st slice drug-ctrl, 2nd slice recov-ctrl
specDiffArea=repmat(nan,[length(bix) n2 2]);
peakFreqDiff=repmat(nan,[length(bix) n2 2]);
sessionName=cell(length(bix),n2);
% in case one or more recordings don't contain one of the analyzed
% behaviors, fill corresponding cell element of mnRdiff with an
% appropriately shaped array of nans (2 slices because currently 2
% difference spectra are calculated)
nanny=repmat(nan,[rfreqix(end),1,2]);
% note that order of loops is different from the one above
for ci=1:n2
  fh=mkfig(ANPAR(1,ci).resPath);
  if size(maIx(ci).rc,1)<20,  orient portrait;
  else  orient tall;
  end
  % temporary var for plotting and derived results calculation
  % indices: freq | channel | conc
  % ** force plotR to have 3 slices even if no recov exists
  plotR=repmat(nan,[rfreqix(end),size(maIx(ci).rc,1),max(n1,3)]);
  Rdiff=repmat(plotR(:,:,1),[1 1 2]);
  % with this definition, all cross measures with princ chan are princ too
  plot_princChInd=find(maIx(ci).rc(:,2)==hp.princChInd);
  for bi=1:length(bix)
    % tell what we're dealing with
    disp(['exprmnt: ' ANPAR(1,ci).resPath ', behavior: ' behav{bi}]);
    % check here, just looking at first channel combo, whether data from one behavior are missing
    if any(any(isnan(R(bi).r{maIx(ci).rc(1,1),maIx(ci).rc(1,2)}(:,:,maIx(ci).sl(1))),1))
      disp(['*** missing data for current behavior, skipping plot']);
      % fill with nans
      for rci=1:size(maIx(ci).rc,1)
        % concatenate nans
        mnRdiff(bi).r{maIx(ci).rc(rci,1),maIx(ci).rc(rci,2)}=...
          cat(2,mnRdiff(bi).r{maIx(ci).rc(rci,1),maIx(ci).rc(rci,2)},nanny);
      end
    else
      for rci=1:size(maIx(ci).rc,1)
        % put in temporary variable for plotting (swapping columns for slices)..
        plotR(:,rci,1:n1)=...
          permute(R(bi).r{maIx(ci).rc(rci,1),maIx(ci).rc(rci,2)}(rfreqix,:,maIx(ci).sl(rci)),[1 3 2]);
      end % for:size(maIx(ci).rc,1)
      % --- plot
      % 1. choose subplot
      subplot(1,nPCols,(bi-1)*2+1);
      % enhance visibility
      set(gca,'color',[.98 .98 .99]);
      % set(gca,'xscale','log','xlim',rfreq([1 end]));
      set(gca,'xlim',rfreq([1 end]));
      % leave some space for title
      rexy('ax',gca,'xfac',1.05,'yfac',.95);
      hold on
      % 2. difference spectra & normalization 
      switch normType
        case 'none'
          % 1st slice: drug-ctrl
          Rdiff(:,:,1)=plotR(:,:,2)-plotR(:,:,1);
          % 2nd slice: wash-ctrl
          Rdiff(:,:,2)=plotR(:,:,3)-plotR(:,:,1);
        case 'ctrlPow'
          % normalize area of control to 1, scale other concs by same number, compute differences
          sar=sum(plotR(normrfreqix,:,1),1)*fResol;
          tmpPlotR=plotR./repmat(sar,[rfreqix(end),1,max(n1,3)]);
          % 1st slice: drug-ctrl
          Rdiff(:,:,1)=tmpPlotR(:,:,2)-tmpPlotR(:,:,1);
          % 2nd slice: wash-ctrl
          Rdiff(:,:,2)=tmpPlotR(:,:,3)-tmpPlotR(:,:,1);
        case 'diffPow'
          % compute differences, normalize area (drug-control) to 1
          % 1st slice: drug-ctrl
          Rdiff(:,:,1)=plotR(:,:,2)-plotR(:,:,1);
          % 2nd slice: wash-ctrl
          Rdiff(:,:,2)=plotR(:,:,3)-plotR(:,:,1);
          sar=sum(abs(Rdiff(normrfreqix,:,1)),1)*fResol;
          Rdiff=Rdiff./repmat(sar,[rfreqix(end),1,2]);
        otherwise
          error('illegal choice of normType');
      end
      
      % analysis
      % note: anfreqix is a local index into rfreqix, which is why Rdiff can be
      % indexed the way it is done below
      % a) integral of difference spectrum 
      tmpsum=sum(Rdiff(anfreqix,plot_princChInd,:),1)*fResol;
      specDiffArea(bi,ci,:)=tmpsum;
      
      % b) frequency of peak
      tmpEv=evdeal(plotR(anfreqix,plot_princChInd,:),'idx',{'minmaxpeak'});
      % prevent indexing errors in case recovery does not exist
      if all(isfinite(tmpEv.maxPeakT([1 2])))      
        peakFreqDiff(bi,ci,1)=diff(rfreq(anfreqix(tmpEv.maxPeakT([1 2]))));
      end
      if all(isfinite(tmpEv.maxPeakT([1 3])))     
        peakFreqDiff(bi,ci,2)=diff(rfreq(anfreqix(tmpEv.maxPeakT([1 3]))));
      end
      sessionName{bi,ci}=[ANPAR(1,ci).resPath ', ' behav{bi} '     '];
        
      % 3. log?
      if logp
        plotR=log10(plotR);
        % log of differences does not make sense
      end

      % 4. compute offsets
      offs=cumsum([0 min(min(plotR(:,1:end-1,:),[],1),[],3)-max(max(plotR(:,2:end,:),[],1),[],3)]);      
      rdoffs=cumsum([0 min(min(Rdiff(:,1:end-1,:),[],1),[],3)-max(max(Rdiff(:,2:end,:),[],1),[],3)]);      
      % 5. plot spectra
      for ri=1:n1
        ph=plot(rfreq,plotR(:,:,ri)+repmat(offs,rfreqix(end),1),'-');
        niceyax; 
        % standard line width
        lwi=1;
        switch ri
          case 1
            app=[];
          case 2
            app='_drug';
            lwi=2.5;
          case 3
            app='_recov';
        end
        % gbix is index into segTypeGlobP
        gbix=strmatch([behav{bi} app],segTypeGlobP(:,1),'exact');
        set(ph,'color',segTypeGlobP{gbix,2},'linewidth',lwi);
        if ri==2
          % now mark principal channels: overplot with thin red line
          ph=plot(rfreq,plotR(:,plot_princChInd,ri)+repmat(offs(plot_princChInd),rfreqix(end),1),'y-');
          set(ph,'linewidth',.5);
        end
      end
      % vertical lines
      line([vlines;vlines],(get(gca,'ylim'))'*ones(size(vlines)),'color',[.7 .7 .7],'linestyle',':');
      title([behav{bi} ', [ctrl,applic,(recov)]']);
      % 6. plot differences of normalized spectra
      subplot(1,nPCols,(bi-1)*2+2);
      % enhance visibility
      set(gca,'color',[.98 .98 .99]);
      set(gca,'xscale','log','xlim',rfreq([1 end]));
      % leave some space for title
      rexy('ax',gca,'xfac',1.05,'yfac',.95);
      hold on
      title('difference spectra');
      % drug-ctrl
      ph=plot(rfreq,Rdiff(:,:,1)+repmat(rdoffs,rfreqix(end),1),'b-');
      % recov-ctrl
      ph=plot(rfreq,Rdiff(:,:,2)+repmat(rdoffs,rfreqix(end),1),'c-');
      niceyax
      % now mark principal channels: overplot with thick line
      ph=plot(rfreq,Rdiff(:,plot_princChInd,1)+repmat(rdoffs(plot_princChInd),rfreqix(end),1),'b-');
      set(ph,'linewidth',2);
      ph=plot(rfreq,Rdiff(:,plot_princChInd,2)+repmat(rdoffs(plot_princChInd),rfreqix(end),1),'c-');
      set(ph,'linewidth',2);
      ultext('blue:drug-ctrl; cyan:recov-ctrl');
      % zero lines (avoid zero as lower xbound in log coordinates)
      line((max(get(gca,'xlim'),fResol))'*ones(1,size(maIx(ci).rc,1)),...
        [rdoffs;rdoffs],'color',[.7 .7 .7],'linestyle','-.')
      % vertical lines
      line([vlines;vlines],(get(gca,'ylim'))'*ones(size(vlines)),'color',[.7 .7 .7],'linestyle',':');
      % 7. use this lazy slot to transfer difference spectra to mnRdiff
      for rci=1:size(maIx(ci).rc,1)
        % concatenate so we have one column per experiment (means & std will be
        % computed further below) 
        mnRdiff(bi).r{maIx(ci).rc(rci,1),maIx(ci).rc(rci,2)}=...
          cat(2,mnRdiff(bi).r{maIx(ci).rc(rci,1),maIx(ci).rc(rci,2)},Rdiff(:,rci,:));
      end 
    end % if:any(isnan(...
  end % for:behav
  % the output of this function is aesthetically extremely repelling, but use it until better stuff
  % available
  figtitle([ANPAR(1,ci).resPath ', norm=' normType]);
  % if respath does not contain a drive letter, pre-pend WP.rootPath 
  if isempty(strfind(ANPAR(1,ci).resPath,':')), ANPAR(1,ci).resPath=[WP.rootPath ANPAR(1,ci).resPath]; end
  if ~isempty(printas)
    print(printas,[ANPAR(1,ci).resPath '\rr_' rv{1} '_' DSET(1,ci).abfFn(1:end-3) '_'...
        int2str(round(rfreq(1))) '_' int2str(round(rfreq(end))) '.jpg']);
  end
end % for:cols of ANPAR=experiments
clear Rdiff plotR tmpPlotR;



% -------- PART VI:  final effort: average differences & plot
% find all common channel combos
mn_rc=unique(cat(1,maIx(:).rc),'rows');
fh=mkfig('means');
if size(mn_rc,1)<20, orient portrait;
else orient tall;
end
% temporary var for plotting, 4D:
% freq | channel | (mean, std) | (drug-ctrl, recov-ctrl)
Rdiff2=zeros(rfreqix(end),size(mn_rc,1),2,2);
% with this definition, all cross measures with princ chan are princ too
plot_princChInd=find(mn_rc(:,2)==hp.princChInd);
% reminder: each cell of mnRdiff shall have this structure (after the loop
% below)
% - freq down columns
% - columns order:  mean | std 
% - slice order: drug-control | recov-control
for bi=1:length(bix)
  Rdiff2(:)=0;
  for rci=1:size(mn_rc,1)
    % compute means and std & additionally put those in plot var
    % (there may be bad columns due to missing behavior)
    if n1==2
      goodColumns=all(isfinite(mnRdiff(bi).r{mn_rc(rci,1),mn_rc(rci,2)}),1);      
    else
      goodColumns=all(all(isfinite(mnRdiff(bi).r{mn_rc(rci,1),mn_rc(rci,2)}),1),3);
    end
    mnRdiff(bi).r{mn_rc(rci,1),mn_rc(rci,2)}=...
      [mean(mnRdiff(bi).r{mn_rc(rci,1),mn_rc(rci,2)}(:,goodColumns,:),2)...
      std(mnRdiff(bi).r{mn_rc(rci,1),mn_rc(rci,2)}(:,goodColumns,:),0,2)];
    % there may be cases where the top- or bottommost channel is present
    % only in one experiment and only one of the behaviors. In this case 
    % goodColumns is logical scalar 0, the corresponding mnRdiff column is
    % all nans and plots don't come out because cumsum messes up. Prevent
    % this by leaving corresponding Rdiff2 columns zero.
    if any(goodColumns)
      Rdiff2(:,rci,:,1)=permute(mnRdiff(bi).r{mn_rc(rci,1),mn_rc(rci,2)}(:,:,1),[1 3 2]);
      Rdiff2(:,rci,:,2)=permute(mnRdiff(bi).r{mn_rc(rci,1),mn_rc(rci,2)}(:,:,2),[1 3 2]);
    end
  end
  % now treat drug-control and recov-drug separately
  for g=1:n1-1
    % compute offsets
    rdoffs=[0 min(Rdiff2(:,1:end-1,1,g)-Rdiff2(:,1:end-1,2,g))-...
      max(Rdiff2(:,2:end,1,g)+Rdiff2(:,2:end,2,g))];
    % replace all zero columns by the average of all 
    rdoffs(~(rdoffs))=mean(rdoffs);
    rdoffs=cumsum(rdoffs);
    % 6. plot differences of normalized spectra
    subplot(1,nPCols,(bi-1)*2+g);
    % enhance visibility
    set(gca,'color',[.98 .98 .99]);
    set(gca,'xscale','log','xlim',rfreq([1 end]));
    hold on
    if g==1, title('drug-control');
    else title('recov-control');
    end
    % plot +/- std first so that they don't overplot means in case of n=1
    ph=plot(rfreq,[Rdiff2(:,:,1,g)+Rdiff2(:,:,2,g) Rdiff2(:,:,1,g)-Rdiff2(:,:,2,g)]+repmat(rdoffs,rfreqix(end),2),'-');
    set(ph,'color',[.7 .7 .7]);
    ph=plot(rfreq,Rdiff2(:,:,1,g)+repmat(rdoffs,rfreqix(end),1),'b-');
    niceyax
    % now mark principal channels: overplot with thick line
    ph=plot(rfreq,Rdiff2(:,plot_princChInd,1,g)+repmat(rdoffs(plot_princChInd),rfreqix(end),1),'b-');
    set(ph,'linewidth',2);
    % zero lines (avoid zero as lower xbound in log coordinates)
    line((max(get(gca,'xlim'),fResol))'*ones(1,size(mn_rc,1)),...
      [rdoffs;rdoffs],'color',[.7 .7 .7],'linestyle','-.')
    % vertical lines
    line([vlines;vlines],(get(gca,'ylim'))'*ones(size(vlines)),'color',[.7 .7 .7],'linestyle',':');
  end
end


figure, hold on
for bi=1:length(bix)
  disp('sessions:')
  disp(sessionName(bi,:)')
  disp('Area of difference spectrum: drug-ctrl | recov-ctrl:')
  permute(specDiffArea(bi,:,:),[2 3 1])
  disp('Peak freq difference (Hz): drug-ctrl | recov-ctrl:')
  permute(peakFreqDiff(bi,:,:),[2 3 1])
  gbix=strmatch(behav{bi},segTypeGlobP(:,1),'exact');
  ph=plot(specDiffArea(bi,:,1),peakFreqDiff(bi,:,1),'o');
  set(ph,'color',segTypeGlobP{gbix,2},'markerfacecolor',segTypeGlobP{gbix,2});
end
nicexyax
grid on
xlabel('Area of difference spectrum (drug-ctrl)');
ylabel('Peak freq difference (drug-ctrl), Hz');

% ----- local funcs -----------

function figha=mkfig(ftag)
figha=findobj('tag',ftag);
if isempty(figha), figha=figure;
else  figure(figha);
end
tmpScrSz=get(0,'Screensize')*.8+.1*rand;
tmpScrSz([1 2])=tmpScrSz([1 2])+tmpScrSz([3 4])*.10;
set(figha,'position',tmpScrSz,'tag',ftag,'name',ftag,...
  'color',[.9 .9 1],'numbertitle','off');
clf;
