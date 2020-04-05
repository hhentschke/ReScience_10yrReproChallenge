function rdeal(rfn,idv,job,varargin)
% ** function rdeal(rfn,idv,job,varargin)
%
%                         >>> INPUT VARIABLES >>>
%
% NAME           TYPE/DEFAULT          DESCRIPTION
% rfn            char arr              name of matfile (without extension)
%                                       containing R, the output of combine_r 
% idv            struct array          contains names and levels of the
%                                       factors in R (independent variables
%                                       e.g. behavior, genotype, etc.) and
%                                       details as to how these shall be
%                                       dealt with in statistical analysis
% job            cell array            job(s) to be performed on the data.
%                                       Currently, 
%   'rm_anova2'       - 2-way ANOVA with repeated measures (both factors r.m.)
%   'anovan'          - n-way ANOVA
%   'curveFit_FTest'  - F-test of fitted curves
% 
% curRv          cell array, {'all'}   results variables to be treated
% plotStyle      char array, 'barhh'   'subjects','barhh' or 'none'
% figSubName     char array, ''        string to be inserted into file name
%                                       of exported graphics
% projSubDir     char array, ''        project subdirectory 
% printas        char array,[]         format for graphics export (any legal
%                                       specifier for matlab 'print'
%                                       function like '-djpeg90'; [] means
%                                       do not export)



% TO DO
% - in cross parameters exclude auto values
% - post-hoc, different 'jobs' (stat analyses, maybe including curve fitting)
% - plots, working style 
% - update documentation (also in collect_*)


% -------------------------------------------------------------------------
%                   I. PRELIMINARIES
% -------------------------------------------------------------------------
curRv={'all'};
plotStyle='barhh';
figSubName='';
projSubDir='';
printas=[];
pvpmod(varargin);

export=0;

curFigPath=['\' projSubDir '\figures'];
close all;
exportFn=['\' projSubDir '\' mfilename '_export.txt'];

% load data 
rmouse_ini;
load([WP.rootPath '\' projSubDir '\' rfn '.mat'],'R','rv','RInfo','Rraw');

% deal with results variable(s) of interest 
if isequal(curRv,{'all'})
  curRv=rv;
  rvix=1:length(rv);
else
  [nada,rvix]=intersect(rv,curRv);
  rvix=sort(rvix);
  if ~isempty(setdiff(curRv,rv))
    warndlg('some rv illegal (forgot to append ''_auto or'' ''_cross''?');
  end
  if isempty(rvix)
    error('no matching rv (forgot to append ''_auto or'' ''_cross''?');
  end
end
% make sure rvix is a row array
if size(rvix,1)>size(rvix,2)
  rvix=rvix';
end

if ~ismember(plotStyle,{'subjects','barhh','none'})
  error('illegal choice of ''plotStyle''');
end


% -------------------------------------------------------------------------
%                   II. SETUP OF FACTORS AND LEVELS
% -------------------------------------------------------------------------
disp('** checking factors and their levels..');
nIdv=length(idv);

% some shorties
recIx=strmatch('rec site',{idv.name});
behIx=strmatch('behavior',{idv.name});
druIx=strmatch('drug',{idv.name});
genIx=strmatch('genotype',{idv.name});
subIx=strmatch('subject',{idv.name});

% *************************************************************************
% although most of the code below does not rely on any particular order of
% factors, the code producing the plots does. Furthermore, there is no
% disadvantage at all of sticking to a particular order. Therefore, check
% that order rigorously
if ~isequal([recIx behIx druIx genIx subIx],1:nIdv),
  error('input structure idv is not set up properly (order and/or number of factors)');
end
% *************************************************************************

% combine_r collects data according to the information contained in DSET,
% ANPAR and RInfo into cell array R.d (and Rraw.d). Here, in rdeal, we may 
% inadvertently want to analyze (as defined by idv) a factor's level not
% present in the compiled data. The lines below check this for genotype,
% behavior and drug. For recording site and subject a different logic
% applies.
checkLevelFac={'behavior','genotype','drug'};
for g=1:length(checkLevelFac)
  tmpLev=setdiff(...
    idv(strmatch(checkLevelFac(g),{idv.name})).level,...
    RInfo(strmatch(checkLevelFac(g),{RInfo.name})).level);
  if ~isempty(tmpLev)
    error(['factor ' checkLevelFac{g} ': level requested for analysis which does not exist in data']);
  end
end

% --- see which factor(s) in the data shall be factor(s) in the statistical
% analysis and for which factor's levels the statistical analysis should be
% carried out separately (sequentially)

% this is the index to the factor(s) to make it into statistical analysis -
% it will be used later on to access columns of the collected data array, 
% coDat
masterAnFacColIx=[];
% this is the index to the factor(s) whose levels shall be analyzed
% separately - it will also be used later for direct access to coDat
masterSepFacColIx=[];

for g=1:nIdv
  % in the code below new fields 'separateVal', 'deleteVal' etc. are
  % generated. they contain the numeric values of the current factor's 
  % level which shall be treated separately ('separateVal'), deleted
  % ('deleteVal'), etc.
  switch idv(g).name
    % rec site needs a special treatment
    case 'rec site'
      % convert cell array to array 
      idv(g).level=cat(1,idv(g).level{:});
      if ~isequal(sort(idv(g).level),idv(g).level)
        error('rec sites must be sorted (0=principal channel, values decreasing in dorsal direction');
      end
      % set up variable depthLim      
      depthLim=idv(g).level([1 end]);
      idv(g).statsVal=(1:length(idv(g).level))';
      idv(g).separateVal=[];
      idv(g).deleteVal=[];
      idv(g).ignoreVal=[];
    otherwise
      idv(g).statsVal=strmatch('stats',idv(g).treat);
      % if only one level of current factor is set to 'stats' reset
      % .statsVal to [] so this factor will be 'ignored' (the data will be
      % used but levels cannot be compared to each other because there is
      % only one)
      if length(idv(g).statsVal)==1
        idv(g).statsVal=[];
      end
      idv(g).separateVal=strmatch('separate',idv(g).treat);
      idv(g).deleteVal=strmatch('delete',idv(g).treat);
      idv(g).ignoreVal=strmatch('ignore',idv(g).treat);
  end
  idv(g).nLevel=length(idv(g).level);
  % check for illegal level properties ('stats' has to be determined again)
  if length([strmatch('stats',idv(g).treat); idv(g).separateVal;...
      idv(g).deleteVal; idv(g).ignoreVal]) ~= idv(g).nLevel
    error([idv(g).name ': at least one level property is illegal (typo?) - must be ''stats'', ''separate'', ''delete'' or ''ignore''']);
  end
  % check for combinations of level properties which are generally illegal 
  % 'stats' and 'separate' cannot be combined
  if length(idv(g).separateVal)>0 && length(idv(g).statsVal)>0
    error([idv(g).name ': ''separate'' and ''stats'' level properties must not be combined']);
  end
  % either all or no level at all must be 'ignore' 
  if length(idv(g).ignoreVal)>0 && length(idv(g).ignoreVal)<idv(g).nLevel 
    error([idv(g).name ': ''ignore'' level property must not be combined with any other property']);
  end
  % ensure that 'subject' is a factor all levels of which are 'ignored'
  % (this ensures proper rearrangement of data for the plots)
  if strmatch(idv(g).name,'subject')
    if length(idv(g).treat) ~= length(idv(g).ignoreVal)
      error('currently, all ''subject'' levels must be set to ''ignore''');
    end
  end
  
  % determine the crucial factors
  if ~isempty(idv(g).statsVal)
    masterAnFacColIx=[masterAnFacColIx g];
  end
  if ~isempty(idv(g).separateVal)
    masterSepFacColIx=[masterSepFacColIx g];
  end
end % for g=1:nIdv
nMasterAnFacColIx=length(masterAnFacColIx);

% check how many factors' levels shall be dealt with in separate analyses
if length(masterSepFacColIx)>1
  error('only one factor''s levels may be treated separately')
end

% if any([idv(masterAnFacColIx).nLevel]==1)
%   error(['within one of the factors only one level is set to ''stats'' property - this defies the purpose of statistical testing']);
% end

% check job
switch job
  case {'rm_anova2'}
    if nMasterAnFacColIx~=2,
      error(['2-way ANOVA needs exactly two factors, currently ' int2str(nMasterAnFacColIx) ' are specified']);
    end
  case {'anovan'}
    if nMasterAnFacColIx<2,
      error(['ANOVAN needs at least two factors, currently ' int2str(nMasterAnFacColIx) ' is specified']);
    end
  case {'curveFit_FTest'}
    if nMasterAnFacColIx~=2,
      error(['F-test of fitted curves needs exactly two factors, currently ' int2str(nMasterAnFacColIx) ' are specified']);
    end
  otherwise
    error('illegal job');
end

% in case plotStyle is 'subjects' find out about the best subplot
% arrangement: assume orient tall and n by n subplots (=ordinate longer
% than abscissa)
if strcmpi(plotStyle,'subjects')
  nPlotCol=floor(sqrt(idv(subIx).nLevel));
  nPlotRow=ceil(idv(subIx).nLevel/nPlotCol);
end
  

% put out some prose as to what will be done
disp('-----------------------------------');
disp(['--- analysis requested: ' job])
disp('--- factors & levels for statistical analysis: ')
for g=masterAnFacColIx
  disp(['- ' idv(g).name ':']);
  if isnumeric(idv(g).level)
    % only rec site is numeric
    disp([repmat('    ',size(idv(g).level,1),1)  num2str(idv(g).level(idv(g).statsVal))]);
  else
    disp([repmat('    ',length(idv(g).statsVal),1)  strvcat(idv(g).level(idv(g).statsVal))]);
  end
end
if ~isempty(masterSepFacColIx)
  disp('--- separate analyses for the following subpopulations: ')
  for g=masterSepFacColIx
    disp(['- ' idv(g).name ':']);
    if isnumeric(idv(g).level)
      disp([repmat('    ',size(idv(g).level,1),1)  num2str(idv(g).level)]);
    else
      disp([repmat('    ',length(idv(g).level),1)  strvcat(idv(g).level)]);
    end
  end
end
  
nDrug=idv(druIx).nLevel;
dvIx=nIdv+1;

% the very outermost loop: results variable
for rvi=rvix
  disp(['******** ' rv{rvi} ' ********']);
  % -------------------------------------------------------------------------
  %               III.A COLLECTION OF DATA & SETUP OF GROUPING VAR
  % -------------------------------------------------------------------------
  % array for collected data and grouping values
  coDat=[];
  % we need as many loops as [potential factors-3] because data for the
  % different rec sites (first factor), drug and subject can be dealt with 
  % without a loop. outermost loop: genotype
  for gi=1:idv(genIx).nLevel
    % next: behavior
    for bi=1:idv(behIx).nLevel
      tmpD=Rraw.d{bi,rvi,gi};
      % number of data values (dep var) in current cell
      tmpN=size(tmpD,1);
      % set up coDat IN ORDER DEFINED BY COMPFAC & append dep var as last
      % column
      % preallocation
      tmpC=repmat(nan,[tmpN*nDrug  nIdv+1]);
      % ** rec site is currently the only factor whose levels are numeric
      % and are therefore used directly instead of grouping values
      tmpC(:,recIx)=repmat(tmpD(:,1),nDrug,1);
      % drug
      tmpC(:,druIx)=reshape(((1:nDrug)'*ones(1,tmpN))',[tmpN*nDrug 1]);
      % behavior
      tmpC(:,behIx)=repmat(bi,tmpN*nDrug,1);
      % genotype
      tmpC(:,genIx)=repmat(gi,tmpN*nDrug,1);
      % subjects
      tmpC(:,subIx)=repmat(Rraw.indv{bi,rvi,gi},nDrug,1);
      % finally, dep var
      tmpC(:,dvIx)=reshape(tmpD(:,2:end),[tmpN*nDrug 1]);
      % append
      coDat=[coDat; tmpC];
    end % for:behav loop
  end % for:genotype loop
  
  % -------------------------------------------------------------------------
  %               III.B TRIMMING OF DATA, ANALYSIS & PLOTS
  % -------------------------------------------------------------------------

  % a logical array pointing to the rows to be kept
  masterRowIx=repmat(logical(1),size(coDat,1),1);
  % --- first loop: delete undesired observations
  for g=1:nIdv
    % special treatment for rec site 
    if g==recIx
      masterRowIx(~ismember(coDat(:,recIx),idv(g).level))=0;
    else
      masterRowIx(ismember(coDat(:,g),idv(g).deleteVal))=0;
    end
  end % for g=1:nIdv
  % mark for deletion entries with NaN (missing data in Rraw)
  masterRowIx=masterRowIx & isfinite(coDat(:,dvIx));
  
  % cut down 
  coDat=coDat(masterRowIx,:);
  % update masterRowIx (will be used below)
  masterRowIx=masterRowIx(masterRowIx);
  
  if export
    if rvi==rvix(1)
      % initialize
      eCoDat=coDat;
    else
      % append dv columns
      eCoDat=[eCoDat coDat(:,dvIx)];
    end
    if ~isequal(eCoDat(:,[recIx,behIx,druIx,genIx,subIx]),...
        coDat(:,[recIx,behIx,druIx,genIx,subIx]))
      error('collection of export data failed')
    end
    if rvi==rvix(end)
      % invert rec depth so that values increase in dorsal direction
      eCoDat(:,recIx)=eCoDat(:,recIx)*-1;
      % save
      save([WP.rootPath exportFn],'eCoDat','-ascii');
      disp(rv(rvix))
    end
  end
  
  % factors to be analyzed separately
  if isempty(masterSepFacColIx)
    nSepFac=1;
  else
    nSepFac=max(1,length(idv(masterSepFacColIx).separateVal));
  end
  
  if ~strcmpi(plotStyle,'none')
    fh=figure(rvi); clf
  end
  % --- second loop: split data up according to the 'separate' factor levels,
  % if requested, and do the statistics
  for g=1:nSepFac
    if isempty(masterSepFacColIx)
      sepFacRowIx=masterRowIx;
    else
      sepFacRowIx=coDat(:,masterSepFacColIx)==idv(masterSepFacColIx).separateVal(g);
      disp(['']);
      disp(['**** ' idv(masterSepFacColIx).level{idv(masterSepFacColIx).separateVal(g)} ':']);
    end
    sepCoDat=coDat(sepFacRowIx,:);
    if length(unique(sepCoDat(:,subIx)))<=1
      warning('too few subjects for statistical analysis - skipping tests');
    else
      % ------------------------------ do it!
      switch job
        case 'rm_anova2'
          stats = rm_anova2(sepCoDat(:,dvIx),sepCoDat(:,subIx),...
            sepCoDat(:,masterAnFacColIx(1)),sepCoDat(:,masterAnFacColIx(2)),...
            {idv(masterAnFacColIx).name})
          % now, if there's an effect in the factor other than the first 
          % (rec site) or an interaction effect rec site x other 
          % investigate all electrodes in pairwise fashion 
          % § this may change in future, have additional input argument(s)
          % to specify exactly how post-hoc tests shall be carried out
          tmpIx=setdiff(masterAnFacColIx,recIx);
          p=[];
          if stats{strmatch(idv(tmpIx).name,stats(:,1),'exact'),end}<=.05 ||...
              stats{strmatch(['rec site x ' idv(tmpIx).name],stats(:,1),'exact'),end}<=.05
            uRecSite=unique(coDat(:,recIx));
            nRecSite=length(uRecSite);
            for gg=1:nRecSite
              ix=sepCoDat(:,recIx)==uRecSite(gg);
              ix1=sepCoDat(:,masterAnFacColIx(2))==1;
              ix2=sepCoDat(:,masterAnFacColIx(2))==2;
              [h,p(gg)]=ttest(sepCoDat(ix&ix1,dvIx),sepCoDat(ix&ix2,dvIx));
            end
            [uRecSite'; p]
          end
        case 'anovan'
          for gg=1:nMasterAnFacColIx
            tmpGroupVar{gg}=sepCoDat(:,masterAnFacColIx(gg));
          end
          [p,table,stats]=anovan(sepCoDat(:,dvIx),tmpGroupVar,...
            'varnames',{idv(masterAnFacColIx).name},...
            'display','off');
          disp(table)
          
%           p=[];
%           if stats{strmatch('drug',stats(:,1),'exact'),end}<=.05
%             for gg=1:nRecSite
%               ix=ds1(:,1)==uRecSite(gg);
%               [h,p(gg)]=ttest(ds1(ix,2),ds2(ix,2));
%             end
%             p
%           end

            % multcompare(stats,'dimension',[1 2 3],'alpha',0.05,'ctype','lsd')
        case 'curveFit_FTest'
          % --- part 1: split data up according to the second independent
          % factor (=the one under investigation)
          ix1=sepCoDat(:,masterAnFacColIx(2))==1;
          ix2=sepCoDat(:,masterAnFacColIx(2))==2;
          % in keeping with code developed a while ago create new variables
          % ds* (ds for data set)
          ds1=sepCoDat(ix1,[recIx dvIx]);
          ds2=sepCoDat(ix2,[recIx dvIx]);
          ds12=sepCoDat(:,[recIx dvIx]);
          % don't forget to invert scale for rec site
          ds1(:,1)=ds1(:,1)*-1;
          ds2(:,1)=ds2(:,1)*-1;
          ds12(:,1)=ds12(:,1)*-1;
          % --- part 2: depending on parameter transform data and set up model
          tmpRv=rv{rvi}(1:strfind(rv{rvi},'_')-1);
          [ft_,fo_,st_,ds1ix,ds2ix,ds12ix]=curveFit2rmousePar(ds1,ds2,ds12,tmpRv);

          % --- part 3: fit and determine quality of fit for wt and ko
          set(fo_,'Startpoint',st_);
          ds1f=fit(ds1(ds1ix,1),ds1(ds1ix,2),ft_ ,fo_);
          ds2f=fit(ds2(ds2ix,1),ds2(ds2ix,2),ft_ ,fo_);
          ds12f=fit(ds12(ds12ix,1),ds12(ds12ix,2),ft_ ,fo_);

          % --- part 4: create curves representing the fits & plot
          % 1. all data pts & fit
          fitx=ds12(1,1):(ds12(end,1)-ds12(1,1))/200:ds12(end,1);
          % cfit objects will be 'fevaluated' automatically
          ds1fit=ds1f(fitx);
          ds2fit=ds2f(fitx);
          ds12fit=ds12f(fitx);

          % colors & symbols: there are as many groups of data as there are
          % levels of the second factor in plotColumnIx, so use these
          % levels' color specs
          plotColumnIx=masterAnFacColIx;
          pCol1=idv(plotColumnIx(2)).pCol{idv(plotColumnIx(2)).statsVal(1)};
          pSymb1=idv(plotColumnIx(2)).pSymb{idv(plotColumnIx(2)).statsVal(1)};
          pCol2=idv(plotColumnIx(2)).pCol{idv(plotColumnIx(2)).statsVal(2)};
          pSymb2=idv(plotColumnIx(2)).pSymb{idv(plotColumnIx(2)).statsVal(2)};

          % working style plots
          if exist('workFh','var')
            figure(workFh)
          else
            workFh=figure;
          end
          clf, hold on
          % first
          ph=plot(ds1(:,1),ds1(:,2),pSymb1);
          set(ph,'color',pCol1);
          ph=plot(fitx,ds1fit,'-');
          set(ph,'color',pCol1);
          % second
          ph=plot(ds2(:,1),ds2(:,2),pSymb2);
          set(ph,'color',pCol2);
          ph=plot(fitx,ds2fit,'-');
          set(ph,'color',pCol2);
          % combo fit in black
          plot(fitx,ds12fit,'k-');
          nicexyax;
          set(gca,'color','m');
          % --- part 5: F-test
          % statistical test for similarity (order: wt crf, ko crf, combo crf)
          [p,F,radj1,radj2]=curvecomp([ds1(ds1ix,:) ds1f(ds1(ds1ix,1))], ...
            [ds2(ds2ix,:) ds2f(ds2(ds2ix,1))],...
            [ds12(ds12ix,:) ds12f(ds12(ds12ix,1))],...
            length(ds1ix)-length(st_),length(ds2ix)-length(st_));
          urtext(['p=' num2str(p,'%1.3f')],.85,'fontsize',12);
          if p<.05
            disp(['***** H0 (identity of depth-response profiles) rejected, p= ' num2str(p)]);
            %       if p<.01, urtext('**',.9,'fontsize',20);
            %       else urtext('*',.9,'fontsize',20);
            %       end
          else
            disp(['H0 (identity of depth-response profiles) not rejected, p= ' num2str(p)]);
            %       urtext('n.s.');
          end
          % information on goodness of fit
          ultext(['R_{adj}(1):' num2str(radj1,3) '; R_{adj}(2):' num2str(radj2,3)]);
          ds1f
          ds2f
          ds12f
          % focus back on main plotting fig
          if exist('fh','var')
            figure(fh);
          end
        otherwise
          error('illegal job')
      end % switch job
    end % if <sufficient number of subjects>
    % invert rec depth so that values increase in dorsal direction
    sepCoDat(:,recIx)=sepCoDat(:,recIx)*-1;
    
    switch plotStyle
      case 'barhh'
        % rearrange data suitable for averaging and plots, always maintaining
        % same order: [recIx druIx genIx behIx]
        plotColumnIx=[recIx druIx genIx behIx];
        plotColumnIx=plotColumnIx(ismember(plotColumnIx,masterAnFacColIx));
        % don't forget to add subject and dep var columns at end
        [rd,x]=regroup(sepCoDat(:,[plotColumnIx subIx end]));
        % this plot makes sense only if the number of factors (excluding 
        % subjects of course) is 2 or 3 and the last(=third) has two levels
        if ~ismember(nMasterAnFacColIx,[2 3]) || idv(plotColumnIx(end)).nLevel~=2
          warning('number of factors or levels of last factor not appropriate for double horizontal bar plot - skipping plot');
        else
          orient portrait
          % positions of compound axes for separate factors (currently only
          % two)
          compAxisPos=[.05 .2 .4 .6; .55 .2 .4 .6];
          % compute means and std across subjects
          % ** relax: owing to the way function 'regroup' re-organizes the
          % data the line below will produce averages separately for each
          % combination of independent factors. In other words, e.g. wt and 
          % ki will be averaged separately, etc. **
          rdMn=nanmean(rd,nMasterAnFacColIx+1);
          rdStd=nanstd(rd,0,nMasterAnFacColIx+1);
          switch nMasterAnFacColIx
            case 2
              [axh,bh]=errbarhh(x,rdMn(:,1),rdMn(:,2),rdStd(:,1),rdStd(:,2),...
                'pos',compAxisPos(g,:));
              % colors of bars: there is only one group of bars in the
              % right and left plot each, so there is no need to use 
              % different colors. Let's have everything in grey
              set(bh,'facecolor',[.6 .6 .6]);
            case 3
              [axh,bh]=errbarhh(x,rdMn(:,:,1),rdMn(:,:,2),rdStd(:,:,1),rdStd(:,:,2),...
                'pos',compAxisPos(g,:));
              % colors of bars: there are as many groups of bars as there are
              % levels of the second factor in plotColumnIx, so use these
              % levels' color specs
              set(bh(1,:),'facecolor',...
                idv(plotColumnIx(2)).pCol{idv(plotColumnIx(2)).statsVal(1)});
              set(bh(2,:),'facecolor',...
                idv(plotColumnIx(2)).pCol{idv(plotColumnIx(2)).statsVal(2)});
              %          if g==1,
              %             legend(bh(:,1),idv(plotColumnIx(end-1)).level(idv(plotColumnIx(end-1)).statsVal),...
              %               'location','SouthOutside');
              %           end
          end % switch nMasterAnFacColIx
          % title, labels
          set(get(axh(1),'title'),'string',...
            [idv(masterSepFacColIx).level{idv(masterSepFacColIx).separateVal(g)} ', ' rv{rvi}]);
          set(get(axh(1),'xlabel'),'string',...
            idv(plotColumnIx(end)).level{idv(plotColumnIx(end)).statsVal(1)});
          set(get(axh(2),'xlabel'),'string',...
            idv(plotColumnIx(end)).level{idv(plotColumnIx(end)).statsVal(2)});
          set(axh(1),'yaxisloc','left','ytick',x);
        end % if: proper number of factors/levels
      case 'subjects'
        orient tall
        labelscale('fontSz',6,'scaleFac',1.0,'lineW',1.0,'markSz',6);
        % ** watch out: subplots will be plotted incrementally for each separately
        % analyzed factor! **
        % rearrange data suitable for plotting individuals, always
        % maintaining same order: [recIx druIx behIx]. genotype can be
        % ignored because we'll plot individuals (which must be of exactly
        % one genotype)
        plotColumnIx=[recIx druIx behIx];
        plotColumnIx=plotColumnIx(ismember(plotColumnIx,masterAnFacColIx));
        % don't forget to add subject and dep var columns at end
        [rd,x]=regroup(sepCoDat(:,[plotColumnIx subIx end]));

        % because the dimension of rd varies according to the number of 
        % indep factors reorder rd into 2D array and plot blocks of columns
        tmpSz=size(rd);
        rd=reshape(rd,[tmpSz(1) prod(tmpSz(2:end))]);
        % the number of columns belonging to one subject
        tmpBlockSz=prod(tmpSz(2:end-1));
        % loop over subjects
        for sIx=1:idv(subIx).nLevel
          subplot(nPlotRow,nPlotCol,sIx), hold on
          plot(x,rd(:,(1:tmpBlockSz)+tmpBlockSz*(sIx-1)),'o-');
          niceyax;
          xlabel('rec depth')
          title(idv(subIx).level{sIx})
        end
        subpax(fh);
    end % switch plotStyle
  end % for g=1:nSepFac (separate factors)
  

  
  if ~isempty(printas)
    if strfind(printas,'ps')
      ext='.ps';
    elseif strfind(printas,'jpeg'),
      ext='.jpg';
    else ext='';
    end
    print(printas,[WP.rootPath curFigPath '\' figSubName '_' rv{rvi} ext]); % '_' idv(idvIx).name(1:4)
  end
end % for: rvi=rvix


% VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
% VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
% VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
return



% ------ PART IV: curve fitting for the sake of statistics 
for rvi=rvix
  figure(rvi), orient landscape
  % loop either over behaviors (comparison of genotypes) or genotypes
  % (comparison of behaviors)
  for loopi=1:length(loopP)  
    % --- part 1: collect data, restrict to recording sites dorsal of &
    % including SLM and invert x axis so that independent var has values
    % >=0 
    % - data sets:
    switch idv(idvIx).name
      case 'genotype'
        ds1=R.d{loopi,rvi,1};
        ds2=R.d{loopi,rvi,2};
        error('''indv'' var not yet set up');
      case 'behavior'        
        ds1=R.d{1,rvi,loopi};
        ds2=R.d{2,rvi,loopi};
        error('''indv'' var not yet set up');
      case 'drug'        
        ds1=R.d{loopi,rvi,1}(:,[1 2]);
        ds2=R.d{loopi,rvi,1}(:,[1 3]);
        indv1=R.indv{loopi,rvi};
        indv2=R.indv{loopi,rvi};
    end
    % restrict depth range
    tmpix=ds1(:,1)<depthLim(1) | ds1(:,1)>depthLim(2);
    ds1(tmpix,:)=[];
    indv1(tmpix)=[];
    tmpix=ds2(:,1)<depthLim(1) | ds2(:,1)>depthLim(2);
    ds2(tmpix,:)=[];
    indv2(tmpix)=[];

    % invert sign
    ds1(:,1)=ds1(:,1)*-1;
    ds2(:,1)=ds2(:,1)*-1;
    % sort
    [ds1,tmpix]=sortrows(ds1,1);
    indv1=indv1(tmpix);
    [ds2,tmpix]=sortrows(ds2,1);
    indv2=indv2(tmpix);
    % - combination of both (only needed here, not for export)
    ds12=sortrows([ds1; ds2]);
    
    % for plots of average data points (that is, the original grand averages) do the same
    % thing (+shuffle & combine columns so that the resulting var has pos|mean|std cols):
    % data sets:
    switch idv(idvIx).name
      case 'genotype'
        mds1=[R.ue{loopi,rvi,1} R.ga{loopi,rvi,1}(:,1:2)];
        mds2=[R.ue{loopi,rvi,2} R.ga{loopi,rvi,2}(:,1:2)];
        titl=[behav{loopi} ', ' rv{rvi} ', ' rType{rvi}];
      case 'behavior'        
        mds1=[R.ue{1,rvi,loopi} R.ga{1,rvi,loopi}(:,1:2)];
        mds2=[R.ue{2,rvi,loopi} R.ga{2,rvi,loopi}(:,1:2)];
        titl=[gt{loopi} ', ' rv{rvi} ', ' rType{rvi}];
      case 'drug'
        mds1=[R.ue{loopi,rvi,1} R.ga{loopi,rvi,1}(:,1:2,1)];
        mds2=[R.ue{loopi,rvi,1} R.ga{loopi,rvi,1}(:,1:2,2)];
        titl=['wt, comprsn drug, ' behav{loopi} ', ' rv{rvi} ', ' rType{rvi}];
    end
    % restrict depth range
    mds1(mds1(:,1)<depthLim(1) | mds1(:,1)>depthLim(2),:)=[];
    mds2(mds2(:,1)<depthLim(1) | mds2(:,1)>depthLim(2),:)=[];
    % invert sign
    mds1(:,1)=mds1(:,1)*-1;
    mds2(:,1)=mds2(:,1)*-1;

    % this is a good spot to complete idv
    idv(1).level=num2str(mds1(:,1));;
    
    
    % --- part 2: depending on parameter transform data and set up model
    [ft_,fo_,st_,ds1ix,ds2ix,ds12ix]=curveFit2rmousePar(ds1,ds2,ds12,rv{rvi});
    
    % --- part 3: fit and determine quality of fit for wt and ko
    set(fo_,'Startpoint',st_);
    ds1f=fit(ds1(ds1ix,1),ds1(ds1ix,2),ft_ ,fo_);    
    ds2f=fit(ds2(ds2ix,1),ds2(ds2ix,2),ft_ ,fo_);
    ds12f=fit(ds12(ds12ix,1),ds12(ds12ix,2),ft_ ,fo_);

    % --- part 4: create curves representing the fits & plot
    % 1. all data pts & fit
    fitx=ds12(1,1):(ds12(end,1)-ds12(1,1))/200:ds12(end,1);
    % cfit objects will be 'fevaluated' automatically
    ds1fit=ds1f(fitx);
    ds2fit=ds2f(fitx);
    ds12fit=ds12f(fitx);

    subplot(length(loopP),2,(loopi-1)*2+1), hold on
    title(titl);
    % first
    ph=plot(ds1(:,1),ds1(:,2),pset{2,1});
    set(ph,'color',pset{1,1});
    ph=plot(fitx,ds1fit,'-');
    set(ph,'color',pset{1,1});    
    % second
    ph=plot(ds2(:,1),ds2(:,2),pset{2,2});
    set(ph,'color',pset{1,2});
    ph=plot(fitx,ds2fit,'-');
    set(ph,'color',pset{1,2});        
    % combo fit in black
    plot(fitx,ds12fit,'k-');
    nicexyax;
    % 2. averages & fit
    subplot(length(loopP),2,(loopi-1)*2+2), hold on
    ph=errorbar(mds1(:,1),mds1(:,2),mds1(:,3),pset{2,1});
    set(ph,'color',pset{1,1});
    ph=plot(fitx,ds1fit,'-');
    set(ph,'color',pset{1,1});    
    ph=errorbar(mds2(:,1),mds2(:,2),mds2(:,3),pset{2,2});
    set(ph,'color',pset{1,2});
    ph=plot(fitx,ds2fit,'-');
    set(ph,'color',pset{1,2});    
    niceyax;
    % --- part 5: F-test
    % statistical test for similarity (order: wt crf, ko crf, combo crf)
    [p,F,radj1,radj2]=curvecomp([ds1(ds1ix,:) ds1f(ds1(ds1ix,1))], ...
                                [ds2(ds2ix,:) ds2f(ds2(ds2ix,1))],...
                                [ds12(ds12ix,:) ds12f(ds12(ds12ix,1))],...
                                length(ds1ix)-length(st_),length(ds2ix)-length(st_));

    urtext(['p=' num2str(p,'%1.3f')],.85,'fontsize',12);
    if p<.05
      disp(['***** H0 (identity of depth-response profiles) rejected, p= ' num2str(p)]); 
%       if p<.01, urtext('**',.9,'fontsize',20);
%       else urtext('*',.9,'fontsize',20);
%       end
    else
      disp(['H0 (identity of depth-response profiles) not rejected, p= ' num2str(p)]);               
%       urtext('n.s.');
    end
    % information on goodness of fit
    ultext(['R_{adj}(1):' num2str(radj1,3) '; R_{adj}(2):' num2str(radj2,3)]);
    ds1f
    ds2f
    ds12f
  end % for: behaviors
  if ~isempty(printas)
    print(printas,[WP.rootPath curFigPath '\' idv(idvIx).name(1:4) '_' rv{rvi} '_' rType{rvi} '_fit.jpg']); 
    % saveas(gcf,[WP.rootPath curFigPath '\' rv{rvi} '_' rType{rvi} '_fit'],'fig');
  end
end % for: parameters