function R=combine_coh(RInfo)
% ** function R=combine_coh(RInfo)
% A variant of combine_r (see below). Deals only with ** results and
% combines them as matrices.
% 
% Help for combine_r:
% The primary job this function does is to loop through a list of *.mat
% files produced by rmouse.m and to combine the data in them. The data
% files to be combined are listed in ANPAR and DSET, global variables which
% must be present when combine_r is called. ANPAR and DSET are struct
% arrays, concatenations of individual and matching(!) AP and DS. They are
% as a rule generated by simple routines in mfiles like collect_*. Their
% layout determines how the data will be combined: (i) one column per
% animal/experimental session (ii) one row per concentration (for
% pharmacological experiments), (iii) one 'slice' (3rd dimension) per
% genotype 
% Furthermore, new parameters derived from raw results may be computed.


global ANPAR DSET

% name of file to save data in 
rfn='rcomb_wtAtropine_cohTh';


% choose electrode depth profile - physical or functional
dp='phys';'func';

% interpolate missing channels?
ipol=0;

% need a shorty for behavior which is inherently a factor in rmouse data
behav=RInfo(strmatch('behavior',{RInfo.name})).level;

% results variables (rv) to collect 
rv={'rawCohMnTh';'gaeCohMnTh'};
% what type are they
rType=repmat({'cross'},length(rv),1);


rmouse_ini;

% -------- PART I: collection of data
[n1 n2 n3]=size(ANPAR);
[dn1 dn2 dn3]=size(DSET);
if ~isequal([n1 n2 n3],[dn1 dn2 dn3]), error('DSET and ANPAR must have equal dimension'); end

% struct holding collected results: R
tmplt=cell(length(behav),length(rv),n3);
tmplt(:)={[]};
% all of the following fields are 3d cell arrays:
% - row=behavior (in order listed above)
% - col=parameter (in order listed above)
% - slice=genotype (in order listed above)
% - each element of the cell array contains this type of data:
% collected data: 2d arr, 1st col electrode pos, 2nd col value (drug
% exprmnts: 3rd+ cols = values ar var. concentrations)
R.d=tmplt;          
% code for individual animal/session (needed for ANOVA with repeated
% measures and may be handy otherwise)
R.indv=tmplt;       
% 1d cell array; for each of the ue, these are the indices into the
% corresponding R.d
R.ueix=tmplt;       
% 'grand average': 2d arr, holding  mean|std|N
R.ga=tmplt;         


% loop over data sets:
% one experiment per column, concentration down the columns, different types of animals
% (e.g. genotypes) in different slices
for i3=1:n3
  for ci=1:n2
    for ri=1:n1
      AP=ANPAR(ri,ci,i3);
      DS=DSET(ri,ci,i3);
      % if one of the (indispensable) fields of DS isempty, the corresponding struct
      % element does not correspond to a data set
      if ~isempty(DS.rawCh)
        rawCh=rmouse_chan;
        % let's make the reasonable assumption that for all data won in one experiment
        % the number of all LFP channes is invariant
        if ri==1
          % template (single data set): 
          tempo=repmat(nan,[nAllLFPCh 1 nAllLFPCh]);
        end
        % load results var..
        if isempty(strfind(AP.resPath,':')), AP.resPath=[WP.rootPath AP.resPath]; end
        load([AP.resPath '\' AP.resFn],'r');
        % ..find behaviors..
        for bi=1:length(behav)
          bix(bi)=strmatch(behav{bi},{r(:).segmentType});
        end
        % ..and extract
        for bi=1:length(bix)
          for rvi=1:length(rv)
            % tell what we're dealing with
            disp(['data: ' AP.resFn ', behavior: ' behav{bi} ', data:' rv{rvi}]);
            % tmpr is the original data
            eval(['tmpr=r(bix(bi)).' rv{rvi} ';']);
            if bi==1 & rvi==1 & ri==1
              % index into R.d rows
              elix=[1:nAllLFPCh]+size(R.d{bi,rvi,i3},1);
            end
            % first data set: set up electrode position in first column
            % slice
            if ri==1
              if strcmpi(dp,'phys')
                % �
                R.d{bi,rvi,i3}(elix,1,1:nAllLFPCh)=repmat(WP.elx,[1 1 nAllLFPCh]);
              elseif strcmpi(dp,'func')
                R.d{bi,rvi,i3}(elix,1,1:nAllLFPCh)=repmat(WP.felx,[1 1 nAllLFPCh]);
              end
            end
            y=tempo;
            % **********************************************************
            % select extraction method depending on results var chosen
            % **********************************************************
            switch rv{rvi}
              case {'rawCohMnTh','gaeCohMnTh'}
                % place contents of 2D cell array elementwise into singleton
                % 3D array
                for iCol=1:nAllLFPCh
                  for iRow=1:nAllLFPCh
                    tmpEl=tmpr{iRow,iCol};
                    if ~isempty(tmpEl)
                      y(iRow,1,iCol)=tmpEl;
                    end
                  end
                end
              otherwise 
                error('current rv not allowed');
            end % switch
            % concatenate all data sets: 1st col electrode pos, 2nd+ col values for
            % concentration in usual order (currently: control, drug, recovery)
            if ~any(any(isfinite(y)))
              warning(['data to be extracted do not exist']);
              R.d{bi,rvi,i3}(elix,ri+1,1:nAllLFPCh)=tempo;
              R.indv{bi,rvi,i3}(elix,1)=(i3-1)*100+ci;
            else
              % now insert
              R.d{bi,rvi,i3}(elix,ri+1,1:nAllLFPCh)=y;
              % individual animal 'code': all wt=column order; all ko=100+column order
              R.indv{bi,rvi,i3}(elix,1)=(i3-1)*100+ci;
            end % if:isempty(y)
          end % for:par
        end % for:behav
      end % if ~isempty(DS.rawCh)
    end % for:rows of ANPAR=concs
  end % for:cols of ANPAR=experiments
end % for:slices of ANPAR=genotypes
% (r not needed anymore)
clear r

% save the 'raw' R with lots of nans in it (missing values)
Rraw=R;

% ------ PART IV: save
% append 'auto' or 'cross' specifier to rv 
for rvi=1:length(rv)
  rv{rvi}=[rv{rvi} '_' rType{rvi}];
end
save(rfn,'Rraw','ANPAR','DSET','RInfo','rv');
