function R=combine3_m
% ** function R=combine3_m
% combines individual rmouse data sets and computes grand averages
% needs ANPAR and DSET (=concatenation of individual and matching AP and DS)
% difference to combine2: 
% - compares drug treatment with control & recovery (matched data; uses DS.conc)
% - has an output variable R
% **** will mess up completely if row order of AP and DS deviates from the following:
%       control
%       application
%       wash (optional)



global ANPAR DSET AP DS WP

% ------------ frequently modified parameters -------------------------------
% choose behaviors to be compared/plotted (legal value of AP.segmentType)
behav={'immobile','exploring'};

% choose auto- or cross-channel results 
q='cross';
q='auto';

% obtain results variable(s) to collect and average/plot
switch q
  case 'auto'
    rv_extended={...
      'detheCCPeakMn', 'delta & theta envelope, peak crosscorrelation, dimensionless';...
      'detheCCPeakTMn', 'delta & theta envelope, lag at peak crosscorrelation, ms';...
      'thgaeCCPeakMn', 'theta & gamma envelope, peak crosscorrelation, dimensionless';...
      'thgaeCCPeakTMn', 'theta & gamma envelope, lag at peak crosscorrelation, ms';...
      'thgaeCCZScore', 'theta & gamma envelope, Z-score of peak crosscorrelation, dimensionless';...
      'thgaeCCZTestP', 'theta & gamma envelope, p-value of peak crosscorrelation, dimensionless';...
      'rawPPeakMn', 'theta, peaks of individual (segmental) power spectra, averaged, mV^2/Hz';...
      'rawPPeakTMn', 'theta, frequency of peaks of individual (segmental) power spectra, averaged, Hz';...
      'rawPMnPeak', 'theta, peak of averaged power spectra, mV^2/Hz';...
      'rawPMnPeakT', 'theta, frequency peak of averaged power spectra, mV^2/Hz';...
      'rawDePEMn', 'delta, power, mV^2';...
      'rawThPEMn', 'theta, power, mV^2';...
      'rawThNarrowPEMn', 'narrow band of theta (peak freq +/-1 Hz), power, mV^2';...
      'rawBePEMn', 'beta, power, mV^2';...
      'rawGaPEMn', 'gamma, power, mV^2';...
      'rawRiPEMn', 'ripples, power, mV^2';...
      'thCCPosPeakDecayMn', 'theta, height of largest side peak in autocorrelation, dimensionless';...
      'thgaeCCPosPeakDecayMn', 'theta & gamma envelope, height of largest side peak in their correlation, dimensionless';...
      'gaePPeakMn', 'gamma envelope, peaks of individual (segmental) power spectra, averaged, mV^2/Hz';...
      'gaePPeakTMn', 'gamma envelope, frequency of peaks of individual (segmental) power spectra, averaged, Hz';...
      'gaePMnPeak', 'gamma envelope, peak of averaged power spectra, mV^2/Hz';...
      'gaePMnPeakT', 'gamma envelope, frequency peak of averaged power spectra, mV^2/Hz';...
      'gaeThPEMn', 'gamma envelope, power in theta band, mV^2';...
      'gaeThNarrowPEMn', 'gamma envelope, power in narrow theta band, mV^2';...
};
  case 'cross'
    rv_extended={...
      'deCCPeakMn', 'delta, peak crosscorrelation, dimensionless';...
      'deCCPeakTMn', 'delta, lag of peak crosscorrelation, ms';...
      'thCCPeakMn', 'theta, peak crosscorrelation, dimensionless';...
      'thCCPeakTMn', 'theta, lag of peak crosscorrelation, ms';...
      'gaCCPeakMn', 'gamma, peak crosscorrelation, dimensionless';...
      'gaCCPeakTMn', 'gamma, lag of peak crosscorrelation, ms';...
      'thHieCCPeakMn', 'theta hi, peak crosscorrelation, dimensionless';...
      'thHieCCPeakTMn', 'theta hi, lag of peak crosscorrelation, ms';...
      'thLoeCCPeakMn', 'theta lo, peak crosscorrelation, dimensionless';...
      'thLoeCCPeakTMn', 'theta lo, lag of peak crosscorrelation, ms';...
      'gaeCCPeakMn', 'gamma envelope, peak crosscorrelation, dimensionless';...
      'gaeCCPeakTMn', 'gamma envelope, lag of peak crosscorrelation, ms';...
      'rawPPeakMn', 'theta, peaks of individual (segmental) power spectra, averaged, mV^2/Hz';...
      'rawPPeakTMn', 'theta, frequency of peaks of individual (segmental) power spectra, Hz';...
      'rawPMnPeak', 'theta, peak of averaged power spectra, mV^2/Hz';...
      'rawPMnPeakT', 'theta, frequency peak of averaged power spectra, mV^2/Hz';...
      'rawCohMnDe', 'delta, coherence, dimensionless';...
      'rawCohMnTh', 'theta, coherence, dimensionless';...
      'rawCohMnThNarrow', 'narrow band of theta, coherence, dimensionless';...
      'rawCohMnBe', 'beta, coherence, dimensionless';...
      'rawCohMnGa', 'gamma, coherence, dimensionless';...
      'rawCohMnRi', 'ripples, coherence, dimensionless';...
      'rawDePEMn', 'delta, power, mV^2';...
      'rawThPEMn', 'theta, power, mV^2';...
      'rawThNarrowPEMn', 'narrow band of theta (peak freq +/-1 Hz), power, mV^2';...
      'rawBePEMn', 'beta, power, mV^2';...
      'rawGaPEMn', 'gamma, power, mV^2';...
      'rawRiPEMn', 'ripples, power, mV^2';...
      'gaeCohMnTh', 'gamma envelope, coherence in theta range, dimensionless';...
      'gaeCohMnThNarrow', 'gamma envelope, coherence in narrow theta range, dimensionless';...
      'gaeCohMnDe', 'gamma envelope, coherence in delta range, dimensionless';...
      };
end


% ------------ rarely modified parameters -------------------------------
% choose electrode depth profile - physical or functional
dp='phys';'func';
% ..and corresponding symbols and colors for plots
pstr={'mo';'gs'};
% shall results be normalized to control (drug-free period)?
% ** watch out, this applies to the computation of R.d; the original (not normalized) 
% values will not be available if you choose this option
normalize=0;


% ------------- do not change anything below here -----------------------------
% ------------- do not change anything below here -----------------------------
rv=rv_extended(:,1);

% struct holding collected results: R
tmplt=cell(length(behav),length(rv));
tmplt(:)={[]};
% all of the following fields are 2d cell arrays:
% - row=behavior (in order listed above)
% - col=parameter (in order listed above)
% - each element of the cell array contains this type of data:
% collected data: 2d arr
% column order: electrode pos | control | applic | wash | conc
R.d=tmplt;

% collection of indices into DSET and ANPAR for display of file names and
% concentrations
R.idx=tmplt;

% char arr identifying individuals
R.indv=tmplt;       
% 1d cell array; for each of the ue, these are the indices into the
% corresponding R.d
R.ueix=tmplt;       
% 'grand average': 2d arr, holding  ElPos|mean|std|N
R.ga=tmplt;         
% this one is a one-row cell array:
% 1d array, for each parameter holding results (p-values etc) 
% from statistical comparison among behaviors 'immobile' and 'exploring'
% ** currently not needed 
R.bstat=tmplt(1,:); 

close all;
labelscale('fontSz',10,'scaleFac',1.0,'lineW',1.0,'markSz',6); 
rmouse_ini;

% -------- PART I: collection of data
[n1 n2 n3]=size(ANPAR);
[dn1 dn2 dn3]=size(DSET);
if ~isequal([n1 n2 n3],[dn1 dn2 dn3]), error('DSET and ANPAR must have equal dimension'); end
% loop over data sets: 
% one experiment per column, concentration down the columns
for ci=1:n2
  for ri=1:n1
    AP=ANPAR(ri,ci);
    DS=DSET(ri,ci);
    if ~isfield(DS,'conc'), error('field conc missing'); end
    rawCh=rmouse_chan;
    % let's make the reasonable assumption that for all data won in one experiment
    % the number of all LFP channes is invariant
    if ri==1
      % template (single data set)
      tempo=repmat(nan,AP.nAllLFPCh,1);
      % template (whole experiment)
      tempo2=repmat(tempo,1,n1);
    end
    % load results var..
    if isempty(strfind(AP.resPath,':')), AP.resPath=[WP.rootPath AP.resPath]; end    
    load([AP.resPath '\' AP.resFn],'r');
    % ..find behaviors..
    tmpb=cell(1,length(r));
    [tmpb{:}]=deal(r.segmentType);
    for bi=1:length(behav)
      bix(bi)=strmatch(behav{bi},tmpb);
    end
    % ..and extract
    for bi=1:length(bix)
      for rvi=1:length(rv)
        % tell what we're dealing with 
        disp(['data: ' AP.resFn ', behavior: ' behav{bi} ', data:' rv{rvi}]);      
        % tmpr is the original data 
        try
            eval(['tmpr=r(bix(bi)).' rv{rvi} ';']);
        catch
            warning([rv{rvi} ' does not exist'])
            tmpr = [];
        end
        if bi==1 & rvi==1 & ri==1
          % index into R.d rows
          elix=[1:AP.nAllLFPCh]+size(R.d{bi,rvi},1);
        end
        % first data set: set up electrode position in first column
        if ri==1
          if strcmpi(dp,'phys')          
            R.d{bi,rvi}(elix,1)=WP.elx;
          elseif strcmpi(dp,'func')
            R.d{bi,rvi}(elix,1)=WP.felx;            
          end
        end

%         % name of experimental animal/session
%         % R.indv{bi,rvi}(elix,1:80)=repmat(sprintf('%-80s',[DS.dpath '; ' num2str(DS.conc)]),length(elix),1);
%         R.indv{bi,rvi}(elix,1:80)=repmat(sprintf('%-80s',DS.dpath),length(elix),1);
        y=[];
        % select extraction method depending on results var chosen 
        switch rv{rvi}
          case {'thgaeCCPeakMn','thgaeCCPeakTMn','detheCCPeakMn','detheCCPeakTMn','thgaeCCPosPeakDecayMn','thgaeCCZScore','thgaeCCZTestP'}
            if ~isempty(tmpr)
              y=tempo;
              if strcmpi(rv{rvi},'thgaeCCZScore')
                % Z>2.5 as criterion (corresponds to ~p=0.013): relative number of segments
                nk=size(tmpr,1);
                tmpr=sum(tmpr>2.5,1)./nk;
              elseif strcmpi(rv{rvi},'thgaeCCZTestP')
                % p<.05 as criterion: relative number of segments
                nk=size(tmpr,1);
                tmpr=sum(tmpr<.05,1)./nk;
              end
              y(AP.LFPccInd)=tmpr;
            end
          otherwise
            if length(celldiag(tmpr))>1
              if strcmpi(q,'cross')
                % extract CC data for principal channel: non analyzed channels 
                % (nans) will be ignored on the plot
                y=cat(1, tmpr{1:AP.LFPpcInd2,AP.LFPpcInd2});
                y=[y; cat(1, tmpr{AP.LFPpcInd2,AP.LFPpcInd2+1:end})];
              elseif strcmpi(q,'auto')
                y=cat(1, tmpr{AP.dixie});              
              end
            end
        end % switch
        % concatenate all data sets: 1st col electrode pos, 2nd+ col values for
        % concentration in usual order (currently: control, drug, recovery)
        if isempty(y)
          warning(['data to be extracted do not exist']);
          R.d{bi,rvi}(elix,ri+1)=tempo;          
        else
          R.d{bi,rvi}(elix,ri+1)=y;
        end % if:isempty(y)
        % put linear index for current element of DSET/ANPAR in position
        % exactly matching its data counterpart in R.d
        R.idx{bi,rvi}(elix,ri+1)=repmat(sub2ind([n1 n2],ri,ci),length(elix),1);
      end % for:par
    end % for:behav
  end % for:rows of ANPAR=concs
end % for:cols of ANPAR=experiments

% ------ PART II: find ue and get indices 
% (this is done separately in order to have grand averaging and statistics run independently of each other)
for bi=1:length(bix)
  for rvi=1:length(rv)
    disp(['behavior: ' behav{bi} ', data:' rv{rvi}]);
    tmpr=R.d{bi,rvi};
    % kick out all electrodes with any nan
    badix=find(any(~isfinite(tmpr),2));
    if ~isempty(badix), disp([int2str(length(badix)) ' entries are NaNs']); end
    tmpr(badix,:)=[];
    R.d{bi,rvi}(badix,:)=[];    
    R.idx{bi,rvi}(badix,:)=[];    
    % electrode positions
    ue=unique(tmpr(:,1));
    % preallocate R.ga
    % dimensions: #electrodes | #parameters (mean,std,N) | #conc (control,drug,recovery)
    R.ga{bi,rvi}=repmat(nan,[length(ue),3,n1]);
    % .. as well as R.bstat
    R.bstat{bi,rvi}=repmat(nan,length(ue),1);
    ix={};
    for ei=1:length(ue)
      ix{ei}=find(tmpr(:,1)==ue(ei));
    end
    R.ueix{bi,rvi}=ix;    
  end
end

% ------ PART III: normalization and data collection
for bi=1:length(bix)
  for rvi=1:length(rv)
    % compute averages and std (nans had been kicked out before)
    tmpr=R.d{bi,rvi};
    % normalization?
    if normalize
      tmpr(:,2:n1+1)=tmpr(:,2:n1+1)./repmat(tmpr(:,2),1,n1);
      R.d{bi,rvi}=tmpr;
    end
    for ei=1:length(R.ueix{bi,rvi})
      ix=R.ueix{bi,rvi}{ei};
      % columns: mean|std|N; slices: control|drug|recovery
      R.ga{bi,rvi}(ei,1,:)=reshape(mean(tmpr(ix,2:n1+1),1),[1 1 n1]);
      R.ga{bi,rvi}(ei,2,:)=reshape(std(tmpr(ix,2:n1+1),0,1),[1 1 n1]);
      R.ga{bi,rvi}(ei,3,:)=repmat(length(ix),[1 1 n1]);      
    end
  end
end
   

% % ------ PART IV: a very sparse dose-response plot
% if n1>1
% for rvi=1:length(rv)
%   figure(rvi), hold on
%   for bi=1:length(bix)
%     % identify specific electrodes
%     tmpi=find(R.d{bi,rvi}(:,1)==0);
%     % plot of original data points vs conc
%     plot(R.d{bi,rvi}(tmpi,end),R.d{bi,rvi}(tmpi,3),pstr{bi});
%   end
%   nicexyax;
% %   xl=get(gca,'xlim');
% %   yl=get(gca,'ylim');
% %   set(gca,'xlim',[0 xl(2)],'ylim',[0 yl(2)]);
%   grid on
%   title([rv{rvi} ', ' q]);
% end
% else
%   warning('dose-response plot not possible because drug data file(s) not specified in DSET/ANPAR');
% end


% ------ Part V:
% extract selected parameters from results structure R created above
disp(' ');
disp(' ');
for bi=1:length(bix)
  disp(' ');
  disp('************************************************************************************');
  disp(['******************** ' behav{bi} ' *************************************************'])
  for pai=1:size(rv_extended,1)
    % current results var
    rr=R.d{bi,pai};
    % princ channel only
    ix=rr(:,1)==0;
    
    % display session names once for each behavior, but don't
    % differentiate between results vars because they should all be the
    % same (check this, though)

    if pai==1
      iidx=R.idx{bi,pai}(ix,2:end);
      indvArr2=[];
      for gg=1:n1
        indvArr=[];
        for hh=1:size(iidx,1)
          indvArr=strvcat(indvArr,[DSET(iidx(hh,gg)).dpath  '\' DSET(iidx(hh,gg)).abfFn '|conc=' num2str(DSET(iidx(hh,gg)).conc)]);
        end
        indvArr2=[indvArr2 repmat(',  ',size(iidx,1),1) indvArr];
      end
      % delete leading commas
      indvArr2(:,1:3)=[];
      disp('Data files|concentrations underlying the data columns below (current behavior only):')
      disp(indvArr2);
      disp('************************************************************************************');
    end
    disp(['- ' rv_extended{pai,2}])
    
    disp([num2str(rr(ix,[2:end]),repmat('%2.6f ',1,n1))]);
    
%     % up to this point, each cell of R.d has n1+1 columns (rec depth | ctrl | app | ...|)
%     % new column(s): values normalized to control condition
%     rr(:,n1+2:n1+2+n1-2)=rr(:,3:3+n1-2)./repmat(rr(:,2),1,n1-1);
%     disp([num2str(rr(ix,[2:end]),repmat('%2.6f ',1,2*n1-1))]);
       
  end
end