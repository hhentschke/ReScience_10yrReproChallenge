function specExtract
% extracts spectra from combinations of data sets
% needs ANPAR and DSET (=concatenation of individual and matching(!) AP and DS)

global ANPAR DSET AP DS

% choose behaviors to be compared/plotted (legal value of AP.segmentType)
behav={'immobile','exploring'};
% ..and corresponding symbols and colors for plots
pstr={'m-';'k-'};
% choose results variable(s) to collect and average/plot - must be a field name of r
% (will be put in eval)
% - peak values (amplitudes & lags) of theta, gamma, and theta-gammaEnv CC
rv={'rawPMn'};

% struct holding collected results: R
shell=cell(length(behav),length(rv));
shell(:)={[]};
% all of the following fields are 2d cell arrays:
% - row=behavior (in order listed above)
% - col=parameter (in order listed above)
% - each element of the cell array contains this type of data:
R.d=shell;          % collected data: 2d arr, 1st col F, 2nd col value
R.fn=shell;         % file names for export

close all;
labelscale('fontSz',10,'scaleFac',1.0,'lineW',1.5,'markSz',9); 
rmouse_ini;

% -------- PART I: collection of data
[n1 n2 n3]=size(ANPAR);
[dn1 dn2 dn3]=size(DSET);
if ~isequal([n1 n2 n3],[dn1 dn2 dn3]), error('DSET and ANPAR must have equal dimension'); end
% loop over data sets: 
% one experiment per column, concentration down the columns
for ci=1:n2
  frix=[];
  % this should not be necessary, but just to make sure..
  R.d=shell;
  R.fn=shell;
  for ri=1:n1
    AP=ANPAR(ri,ci);
    DS=DSET(ri,ci);
    if ~isfield(DS,'conc'), error('field conc missing'); end
    rmouse_APcheck;
    rawCh=rmouse_chan;
    if ri==2
      R.name{ci}=[DS.dpath '; ' num2str(DS.conc)];
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
        % first data set: set up F in first column and restrict range
        if ri==1
          frix=r(1).F<200;
          R.d{bi,rvi}=r(1).F(frix);
        end
        % first data set: file name for export
        if ri==1
          % kick out all but topmost subdir 
          slix=findstr(DS.dpath,'\');
          R.fn{bi,rvi}=[rv{rvi} '_' DS.dpath(slix(end)+1:end) '_'  behav{bi} '.txt'];
        end
        eval(['y=r(bix(bi)).' rv{rvi} ';']);
        % concatenate data sets: 1st col F, 2nd+ col values for
        % concentration in usual order (currently: control, drug, recovery)
        if length(diag(y))>1        
          R.d{bi,rvi}(:,ri+1)=y{AP.LFPpcInd2,AP.LFPpcInd2}(frix);
        else
          warning(['data to be extracted do not exist']);
          R.d{bi,rvi}(:,ri+1)=nan;
        end % if:length(diag(y))>1
      end % for:par
    end % for:behav
  end % for:rows of ANPAR=concs
  % export to ASCII
  for bi=1:length(bix)
    for rvi=1:length(rv)
      dd=R.d{bi,rvi};
      save(R.fn{bi,rvi},'dd','-ASCII');
    end
  end
end % for:cols of ANPAR=experiments



% ------ PART x: plot
for rvi=1:length(rv)
  figure(rvi), hold on
  for bi=1:length(bix)
    subplot(2,1,bi)
    ph=plot(R.d{bi,rvi}(:,1),R.d{bi,rvi}(:,2:end));
    niceyax;
    % set(gca,'xscale','log','yscale','log');
    set(gca,'xscale','log');
    xl=get(gca,'xlim');
    yl=get(gca,'ylim');
    set(gca,'xlim',[1 100]);
    title([rv{rvi} ', ' behav{bi}]);
  end  
  legend(ph,{'ctrl','drug','recov'})
end
     
