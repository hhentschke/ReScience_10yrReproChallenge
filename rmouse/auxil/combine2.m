function combine2
% combines individual rmouse data sets and computes grand averages
% needs ANPAR and DSET (=concatenation of individual and matching AP and DS)
% - statistics: comparison of behaviors (ttest without correction)

global ANPAR DSET

% choose electrode depth profile - physical or functional
dp='phys';'func';
% choose auto- or cross-channel results (cross not computed for theta, gamma env corr)
q='cross';
q='auto';
% choose behaviors to be compared/plotted (legal value of AP.segmentType)
behav={'immobile','exploring'};
% ..and corresponding symbols and colors for plots
pstr={'mo';'gs'};
% choose results variable(s) to collect and average/plot - must be a field name of r
% (will be put in eval)
switch q
  case 'auto'
    % all that make sense as autos 
    rv={'detheCCPeakMn','detheCCPeakTMn',...
      'thgaeCCPeakMn','thgaeCCPeakTMn',...
      'rawPPeakMn','rawPPeakTMn',...      
      'rawPMnPeak','rawPMnPeakT',...
      'rawDePEMn','rawThPEMn','rawGaPEMn','rawRiPEMn'};
  case 'cross'
    % all that make sense as cross-measures
    rv={'deCCPeakMn','deCCPeakTMn',...
      'thCCPeakMn','thCCPeakTMn',...
      'gaCCPeakMn','gaCCPeakTMn',...
      'thHieCCPeakMn','thHieCCPeakTMn',...      
      'thLoeCCPeakMn','thLoeCCPeakTMn',...      
      'gaeCCPeakMn','gaeCCPeakTMn',...
      'rawPPeakMn','rawPPeakTMn',...      
      'rawPMnPeak','rawPMnPeakT',...
      'rawCohMnDe','rawCohMnTh','rawCohMnGa','rawCohMnRi',...
      'rawDePEMn','rawThPEMn','rawGaPEMn','rawRiPEMn'};
end

% check dimension of DSET and ANPAR
[n1 n2 n3]=size(ANPAR);
[dn1 dn2 dn3]=size(DSET);
if ~isequal([n1 n2 n3],[dn1 dn2 dn3]), error('DSET and ANPAR must have equal dimension'); end
if all([n1 n2]>1)
  error([mfilename ' requires single-column or single-row DSET and ANPAR']);
end

% struct holding collected results: R
R.type='wt';      % what kind of beasts?
tmplt=cell(length(behav),length(rv));
tmplt(:)={[]};
% all of the following fields are 2d cell arrays:
% - row=behavior (in order listed above)
% - col=parameter (in order listed above)
% - each element of the cell array contains this type of data:
R.d=tmplt;          % collected data: 2d arr, 1st col electrode pos, 2nd col value
R.ueix=tmplt;       % 1d cell array; for each of the ue, these are the indices into the corresponding R.d
R.ga=tmplt;         % 'grand average': 2d arr, holding  ElPos|mean|std|N
% this one's a one-rower
R.bstat=tmplt(1,:); % 1d array, for each parameter holding results (p-values etc) 
                    % from statistical comparison among behaviors 'immobile' and 'exploring'

close all;
labelscale('fontSz',10,'scaleFac',1.0,'lineW',1.0,'markSz',6); 

rmouse_ini;

% -------- PART I: collection of data
% loop over data sets
for ii=1:length(ANPAR)
  AP=ANPAR(ii);
  DS=DSET(ii);
  rmouse_APcheck;
  rawCh=rmouse_chan;
  % a template
  tempo=repmat(nan,nAllLFPCh,1);
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
      % select extraction method depending on results var chosen 
      y=[];
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
          if length(diag(tmpr))>1
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
      % concatenate all data sets: 1st col electrode pos, 2nd col value
      if isempty(y)
        warning(['data to be extracted do not exist']);
        y=tempo;
      else
        if strcmpi(dp,'phys')
          R.d{bi,rvi}=cat(1,R.d{bi,rvi},[WP.elx y]);
        elseif strcmpi(dp,'func')
          R.d{bi,rvi}=cat(1,R.d{bi,rvi},[WP.felx y]);
        end
      end % if:isempty(y)
    end % for:par
  end % for:behav
end % for:length(ANPAR)

% ------ PART II: find ue and get indices 
% (this is done separately in order to have grand averaging and statistics run independently of each other)
for bi=1:length(bix)
  for rvi=1:length(rv)
    tmpr=R.d{bi,rvi};
    % find electrodes with any nan
    badix=find(any(~isfinite(tmpr),2));
    if ~isempty(badix),
      disp([int2str(length(badix)) ' entries are NaNs']);
      tmpr(badix,:)=[];
      R.d{bi,rvi}(badix,:)=[];
    end
    % electrode positions
    ue=unique(tmpr(:,1));
    % preallocate R.ga
    r=repmat(nan,length(ue),4);
    % 1st column=electrode positions
    r(:,1)=ue;
    R.ga{bi,rvi}=r;
    % preallocate R.bstat
    if bi==1
      R.bstat{1,rvi}=repmat(nan,length(ue),1);
    end
    ix={};
    for ei=1:length(ue)
      ix{ei}=find(tmpr(:,1)==ue(ei));
    end
    R.ueix{bi,rvi}=ix;    
  end
end

% ------ PART III: grand averages
for bi=1:length(bix)
  for rvi=1:length(rv)
    % compute averages and std:
    tmpr=R.d{bi,rvi};
    % nans had been kicked out before
    for ei=1:length(R.ueix{bi,rvi})
      ix=R.ueix{bi,rvi}{ei};
      % electrode pos|mean|std|N
      R.ga{bi,rvi}(ei,2)=mean(tmpr(ix,2));
      R.ga{bi,rvi}(ei,3)=std(tmpr(ix,2));      
      R.ga{bi,rvi}(ei,4)=length(ix);      
    end
  end
end
   

% ------ PART IV: statistics - two behaviors only
if length(bix)~=2,
  warning('statistics not computed - number of behaviors to compare must be=2');
else
  % first thing to do: check whether each animal has data for the two different behaviors,
  % in which case the data are matched
  for rvi=1:length(rv)  
    if isequal(R.d{1,rvi}(:,1),R.d{2,rvi}(:,1)) & isequal(cat(1,R.ueix{1,rvi}{:}),cat(1,R.ueix{2,rvi}{:}))
      tmpr=[R.d{1,rvi}(:,2) R.d{2,rvi}(:,2)];
      bi=1;
      for ei=1:length(R.ueix{bi,rvi})
        ix=R.ueix{bi,rvi}{ei};
        if length(ix)>1
          [H,p]=ttest(tmpr(ix,1),tmpr(ix,2));
          R.bstat{bi,rvi}(ei)=p;             
        end
      end
    end
  end
end



% ------ PART IV: plot
for rvi=1:length(rv)
  figure(rvi), hold on
  for bi=1:length(bix)
    ph{bi}=errorbar(R.ga{bi,rvi}(:,1),R.ga{bi,rvi}(:,2),zeros(size(R.ga{bi,rvi}(:,2))),R.ga{bi,rvi}(:,3),pstr{bi});
  end
  niceyax;
  % plot on top: crosses for nans, circles for p>=.05, single star for p<.05, double star for p<.01
  yl=get(gca,'ylim');
  yp=yl(2)-diff(yl)*[.05 .1];
  nanix=find(~isfinite(R.bstat{1,rvi}));
  nsix=find(R.bstat{1,rvi}>=.05);
  ix01=find(R.bstat{1,rvi}<.01);  
  ix05=setdiff(find(R.bstat{1,rvi}<.05),ix01);
  plot(R.ga{bi,rvi}(nanix,1), repmat(yp(1),size(nanix)),'kx');
  plot(R.ga{bi,rvi}(nsix,1), repmat(yp(1),size(nsix)),'ko'); 
  plot(R.ga{bi,rvi}(ix05,1), repmat(yp(1),size(ix05)),'k*');   
  plot(R.ga{bi,rvi}(ix01,1), repmat(yp(1),size(ix01)),'k*');     
  plot(R.ga{bi,rvi}(ix01,1), repmat(yp(2),size(ix01)),'k*');     
  % restrain x axis range
  title([rv{rvi} ', ' q]);
end
    



