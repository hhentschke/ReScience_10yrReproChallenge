function varargout=rmouse(varargin)
% ** function varargout=rmouse(varargin)
% Main function for the set of routines analyzing hippocampal field potential
% data from behaving rodents (e.g. running mice, hence the name). See
% documentation for how it works. 
% The input parameters listed below are optional and must be specified as
% parameter/value pairs, e.g. as in 
%          rmouse('callJob','seg_thetaTrigStream');
%                         >>> INPUT VARIABLES >>>
%
% NAME              TYPE/DEFAULT          DESCRIPTION
% masterAP          char arr, []          name of a analysis parameter (AP) file 
%                                         which will be run (eval-ed) before analysis 
%                                         proper starts. This is a convenient way to
%                                         define parameters which should be equal for 
%                                         all data sets to be analyzed
% af_tsl_ext        double array, []      a time stamp list of artifacts
%                                         (determined outside rmouse). It will 
%                                         only be effective if job 'det_artifacts' 
%                                         is specified (all time stamps will be combined).
%                                         ** NOTE: t=0 in af_tsl_ext MUST
%                                         be the beginning of the
%                                         recording, not the beginning of
%                                         the excerpt (AP.rawExcerpt)!!
% callJob           char arr, []          name of a function or script to be
%                                         called (once all other jobs have been
%                                         dealt with). Like masterAP, this must
%                                         be an expression which can be put into
%                                         the eval function. This option is useful 
%                                         for running specialized scripts/functions 
%                                         which rely on the preparatory work done in 
%                                         rmouse but are not part of the
%                                         standard staple of analysis jobs
%                     
%                         <<< OUTPUT VARIABLES <<<
%
% NAME             TYPE/DEFAULT               DESCRIPTION
% varargout{1}     extended time stamp list   time stamps for behavioral
%                                             transitions WITHOUT taking artifacts 
%                                             or too brief segments into account. 
%                                             See etslconst.m for further 
%                                             information on extended time 
%                                             stamp lists.
% varargout{2}     cell array of char arrays  strings describing the behaviors.
%                                             The order in which they are listed
%                                             in this cell array corresponds to
%                                             their numerical code in
%                                             above-mentioned etsl
% varargout{3}     extended time stamp list   final time stamp list for behavioral
%                                             transitions taking artifacts 
%                                             and too brief segments into account. 

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   CHECK OF DATA & PREPARATIONS
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
global DS AP WP D r logstr

masterAP=[];
callJob=[];
af_tsl_ext=[];

pvpmod(varargin);

% --- if provided as input argument, run the AP file
if ~isempty(masterAP)
  if ischar(masterAP)
    try
      eval([masterAP ';']);
    catch
      % if it didn't work, stop
      warning('Breaking because setup of master analysis parameter structure AP failed:');
      disp(lasterr);
      return
    end
  else
    warning(['input parameter to ' mfilename ' must be a char array']);
    return
  end
end

  
% ------ set up a few essential 'working parameters'
rmouse_ini;

% ---- set up logging to file
% check for paths: this must be done before logging because the logfile will reside
% in one of them
% if dpath does not contain a drive letter, pre-pend WP.rootPath
if isempty(strfind(DS.dpath,':')), DS.dpath=[WP.rootPath DS.dpath]; end
% same with results path & stream dir
if isempty(strfind(AP.resPath,':')), AP.resPath=[WP.rootPath AP.resPath]; end
if isempty(strfind(AP.strmDir,':')), AP.strmDir=[WP.rootPath AP.strmDir]; end
if ~exist(DS.dpath,'dir'), error('data path (DS.dpath) does not exist'); end
if ~exist(AP.resPath,'dir'), error('results path (AP.resPath) does not exist'); end
% set up log string and file detailing analysis results and warnings 
logF=[AP.resPath '\' AP.logFn];
diary(logF)
disp(' '); disp(' '); disp(' '); 
% set up the string - and display it so diary puts it in the logfile
logstr={['********** data file: ' DS.abfFn ', analysis started ' datestr(now)]};
disp(logstr{end});

% ---- make sure the current AP has all required fields 
disp('checking fields of AP');
rmouse_apcheck;

% ------ set up channel indices
rawCh=rmouse_chan;

% ------ deal with behavioral issues
% IMPORTANT: sort AP.behavType according to trigger level
% order: [least negative..most negative, least positive..most positive]
mi=[];
tmp1=cat(1,AP.behavType{:,2});
% negative levels
tmpi2=find(tmp1(:,1)<0);
if ~isempty(tmpi2), 
  [nix,tmpi4]=sort(abs(tmp1(tmpi2,1)));
  mi=tmpi2(tmpi4);
end
% positive levels
tmpi3=find(tmp1(:,1)>=0);
if ~isempty(tmpi3), 
  [nix,tmpi5]=sort(tmp1(tmpi3,1));
  mi=[mi; tmpi3(tmpi5)];
end
AP.behavType=AP.behavType(mi,:);
clear tmp* mi nix
% generate new field 'segmentType' which is a generalization of behavType.
% It comprises additional segment types, to be assigned on the basis 
% of neuronal data or timing or whatever
AP.segmentType=cat(1,AP.behavType,{...
    'tooBrief',nan, [.5 .9 .8], '<';...
    'artifact', nan, [1 1 0], 'd';...
    'undef',nan, [.7 .7 .7], '*'...  
});
% more constants needed below
nSegmentType=size(AP.segmentType,1);
% -- create (integer) variables whose value codes for segment type
% -> once again let it be stated that the whole workings of the program depends on
% the (row) order of AP.segmentType and .behavType, which must not be touched
% after this point
immV=strmatch('immobile',AP.segmentType(:,1),'exact');
explV=strmatch('exploring',AP.segmentType(:,1),'exact');
grooV=strmatch('grooming',AP.segmentType(:,1),'exact');
badV=strmatch('bad',AP.segmentType(:,1),'exact');
briefV=strmatch('tooBrief',AP.segmentType(:,1),'exact');
artifV=strmatch('artifact',AP.segmentType(:,1),'exact');
undefV=strmatch('undef',AP.segmentType(:,1),'exact');
disp('*** behavioral codes:');
for g=1:nSegmentType
  disp(sprintf('%20s: %2.0f',AP.segmentType{g,1},g));
end
% -- trigger thresholds of different types of behavior - all in one, sorted in
% ascending order
trigThresh=cat(1,AP.behavType{:,2});

% ------ miscellanea 
% make printas a cell array if it isn't one already
if ~iscell(AP.printas)
  AP.printas={AP.printas};
end

% 'local' variables - saves some typing effort
dpath=DS.dpath;
abfFn=DS.abfFn;
bsFn=AP.bScoreFn;
verbose=WP.verbose;
% (global) variables needed for dealing with extended time stamp lists
etslconst;
% temporary column to hold stop time stamps of segments 
stopCol=etslc.nCol+1;
% some handy constants for placing subplots on summary figures
WP.ymarg=.035;
WP.xmarg=.02;

% display AP & DS once for the logfile
disp(DS);
disp(AP);

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   WORK ON ANALYSES AS SPECIFIED IN AP.job
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
done=0;
% insert initial precomputation job
AP.job(2:end+1)=AP.job(1:end);
AP.job{1}='precomp';
% if callJob is nonempty and a string, insert it
if ~isempty(callJob)
  if ischar(callJob)
    if strcmpi(callJob(1:4),'seg_')
      AP.job{end+1}='seg_callJob';
    else
      AP.job{end+1}='callJob';
    end
  else
    warning('input argument callJob must be a string');
  end
end
% if any of the jobs computing theta phase (lag) information is called, 
% run the routine converting phase lags from ms to rad
[nada,tmpi]=intersect(AP.job,{'seg_rawSpecPow','seg_thetaCC','seg_gammaEnvCC','seg_thetaGammaEnvCC'});
if ~isempty(tmpi),
  AP.job(tmpi(end)+2:end+1)=AP.job(tmpi(end)+1:end);
  AP.job{tmpi(end)+1}='seg_cc_ms2rad';
end
% finally, if any of the jobs computing behaviorally scored segments is
% requested, insert the job preparing readin and preallocation of results
% variables. This also applies to callJob.
tmpi=strmatch('seg_',AP.job);
% however, if job 'segPrep' has been (manually, sorts of) included by user
% don't run it twice
if ~isempty(tmpi) && isempty(strmatch('segPrep',AP.job)),
  AP.job(tmpi(1)+1:end+1)=AP.job(tmpi(1):end);
  AP.job{tmpi(1)}='segPrep';
end

while ~done
  partJob=AP.job{1};
  switch partJob
    case 'precomp'
      % variables which must exist but may be empty
      etsl=[];
      r_etsl=[];
      af_etsl=[];
      % before going about the job, do some checks:
      % AP.bbor
      if ~isfinite(AP.bbor)
        logstr{end+1}='illegal value for AP.bbor';
        error(logstr{end});
      end
      % existence of neuronal data - if matfile exists, pick it instead of abf file
      if exist([dpath '\' abfFn '.mat'],'file')
        load([dpath '\' abfFn '.mat'],'fi');      
      elseif exist([dpath '\' abfFn '.abf'],'file')
        [nix,nix2,fi]=abfload([dpath '\' abfFn '.abf'],'info');
      elseif exist([dpath '\' abfFn '.raw'],'file')
        rawFInfo=rawinfo('filename',[dpath '\' abfFn '.raw'],'print','no');
        fi=rawfi2genfi(rawFInfo,'channels',AP.rawChAnNm);
      else
        error([dpath '\' abfFn ' does not exist'])
      end
      % original sampling interval
      WP.osi=fi.si;      
      % end & beginning of recording in ms
      WP.eor=1e-3*fi.dataPtsPerChan*WP.osi;
      WP.bor=0;
      % end & beginning of excerpt investigated in ms
      if ischar(AP.rawExcerpt)      
        if strcmpi(AP.rawExcerpt,'full'),
          WP.eoe=WP.eor;
          WP.boe=WP.bor;
        else
          error('check AP.rawExcerpt');          
        end
      else 
        WP.eoe=AP.rawExcerpt(2)*1000;
        WP.boe=AP.rawExcerpt(1)*1000;
      end
      % behavioral data: check which type
      fex1=exist([dpath '\' bsFn '.abf'],'file');
      fex2=exist([dpath '\' bsFn '.txt'],'file');
      if fex1 && fex2
        fex1=0;
        logstr{end+1}='found both *.abf and *.txt files for behavioral scoring data - using the latter';
        warning(logstr{end});
      end
      if fex1
        [nix,nix2,b_abfi]=abfload([dpath '\' bsFn '.abf'],'info');
        % b_abfi=abfinfo([dpath '\' bsFn '.abf'],'verbose',verbose);      
        % if length of recordings differ by more than 2 % issue a warning
        tmp1=abs(1-b_abfi.dataPtsPerChan/fi.dataPtsPerChan);
        if tmp1>.02,
          logstr{end+1}='abf files for behavioral scoring and neuronal data differ in length by more than 2 %';
          warning(logstr{end});
        end
        % ms, length of behavioral recording
        b_eor=1e-3*b_abfi.dataPtsPerChan*b_abfi.si;
        % in the time frame set by the neural recording, this is the actual
        % offset to be added to the behav time stamps (ms)
        b_delta=1000*AP.bbor;    
        % start and stop points of excerpt to read from behav file 
        % within its own time frame (ms) 
        bstart=0;
        bstop=min(b_eor,WP.eoe-1000*AP.bbor);
      elseif fex2
        % load 
        bsl=textread([dpath '\' bsFn '.txt'],'','commentstyle','matlab');              
        % the last column should contain zeroes (end of line), so check &
        % delete it
        if size(bsl,2)~=2
          error([dpath '\' bsFn '.txt must contain two columns']);
        end
        % convert from min.sec to ms: 
        bsl(:,1)=floor(bsl(:,1))*60000+(bsl(:,1)-floor(bsl(:,1)))*1e5;
        % neither sort nor do any other thing to list here
      else
        logstr{end+1}='behavioral scoring data: found neither *.abf nor *.txt file)';
        disp(logstr{end});
        error(logstr{end});
      end
      % in title string replace underscores (tex interpreter makes text to subscript)
      tmpfn1=DS.abfFn; tmpfn1(strfind(tmpfn1,'_'))='-'; 
      tmpfn2=AP.bScoreFn; tmpfn2(strfind(tmpfn2,'_'))='-';         
      % title of results figures
      WP.figTitle=[DS.aName ' (' tmpfn1 ' + ' tmpfn2 ')'];
      % core name of results graphics file(s) 
      WP.figName=[AP.resPath '\' AP.resFn];
      % generate the primary analysis results figure
      WP.sumFigH=mkfig('specSum'); 
      orient landscape
      labelscale('fontSz',7,'scaleFac',1.0,'lineW',1.0,'markSz',4); 
      drawnow
      clear tmp* xtl* fex*
      AP.job(1)=[];

    % detect artifacts in neural signals & generate time stamp list    
    case 'det_artifacts'
      if length(AP.rawChMonNm) ~= length(AP.afThresh) || isempty(AP.rawChMonNm)
        logstr{end+1}=['data file: ' AP.resPath '\' DS.abfFn '; scoring file: ' AP.resPath '\' AP.bScoreFn   ' - job ''det_artifacts'' was requested, but thresholds/channels are empty or nonmatching in number'];
        disp(logstr{end});
        warndlg(logstr{end});
      else
        logstr{end+1}='---------- detecting artifacts in neuronal data..';
        disp(logstr{end});
        af_etsl=rmouse_detartifact(artifV,af_tsl_ext);
      end 
      % if next job is not 'gen_btsl', set to internal job 'gen_bseq' 
      jidx=strmatch('gen_btsl',AP.job);
      if isempty(jidx)
        AP.job{1}='gen_bseq';
      elseif jidx==2
        AP.job(1)=[];
      else
        error('the sequence of jobs specified in AP.job is not correct');
      end
      
    % generate time stamp list of behavioral segments  
    case 'gen_btsl'
      if exist('bsl','var') 
        logstr{end+1}='---------- detecting behavioral triggers from txt file..';
        % here, all triggers will be read in, irrespective of values of WP.boe and WP.eoe,
        % and will be cropped further below
        disp(logstr{end});        
        % sort & 3rd column for behavType values
        bsl=sortrows(bsl);
        bsl(:,3)=nan;
        % assign behavType values to trigger levels
        for i=1:size(trigThresh,1)
          triggi=sort(trigThresh(i,:));
          bsl(bsl(:,2)>=triggi(1) & bsl(:,2)<triggi(2),3)=i;
        end
        if any(isnan(bsl(:,3))),
          logstr{end+1}=[dpath '\' bsFn '.txt contains illegal trigger values'];
          disp(logstr{end});
          error(logstr{end});        
        end
        % now copy to etsl
        etsl=bsl(:,1);
        etsl(:,etslc.tagCol)=bsl(:,3);
        % *** set time frame: add AP.bbor and that's it!
        etsl(:,etslc.tsCol)=etsl(:,etslc.tsCol)+1000*AP.bbor;
        clear bsl triggi
      else
        logstr{end+1}='---------- detecting behavioral triggers from ABF file..';
        disp(logstr{end});
        etsl=rmouse_genbetsl(bstart,bstop,b_delta,trigThresh);
      end

      % shift time stamps by AP.bsRt
      etsl(:,etslc.tsCol)=etsl(:,etslc.tsCol)-AP.bsRt;
      % the task of the following lines is to 
      % - delete time stamps outside the range of the requested excerpt
      % - to add a fake time stamp (ts=0) to the etsl's head 
      % because every single tick of raw data must be assigned to one of the segment types.
      % the fake time stamp is either of the same type as the one right before WP.boe 
      % or of type 'undef'.
      % ** etsl must have been sorted at this point **
      killIdx=find(etsl(:,etslc.tsCol)<WP.boe);
      if isempty(killIdx)
        % the 'if find..' case below ensures that appending of an undef time stamp is
        % prevented in case there is a time stamp at exactly 0 (as may be the case in
        % text scoring data)
        if find(etsl(1,etslc.tsCol))
          % in this case, time stamp at t=WP.boe set to undefined (because that is what the
          % time period to the first 'real' behavioral time stamp must be)
          etsl=[zeros(size(etsl(1,:))); etsl];
          etsl(1,[etslc.tsCol etslc.tagCol])=[WP.boe undefV];
        end
      else
        % the last trigger before WP.boe: make the fake trigger of same type as its
        % preceding one
        etsl(killIdx(end),etslc.tsCol)=WP.boe;
      end
      % now delete time stamps outside the interval of neural data recording
      killIdx=[killIdx(1:end-1); find(etsl(:,etslc.tsCol)>=WP.eoe)];
      etsl(killIdx,:)=[];
      logstr{end+1}=[int2str(length(killIdx)) ' triggers beyond extent of neural data currently analyzed were cleared from behavioral event list'];
      disp(logstr{end});
      % deal with ambiguous/erroneous data ('r' = 'raw' in the sense of not corrected):
      r_etsl=etsl;
      % 1. two or more consecutive triggers of equal impulse level
      tmp=~(diff(etsl(:,etslc.tagCol)));
      etsl([false; tmp],:)=[];
      logstr{end+1}=[int2str(length(find(tmp))) ' redundant triggers cleared from behavioral event list'];
      disp(logstr{end});
      % fill the 'duration' column 
      etsl(:,etslc.durCol)=diff([etsl(:,etslc.tsCol); WP.eoe]);      
      % (keep a copy of etsl at this point for continuous time course plots of e.g. 
      % theta peaks. don't confuse with petsl. also, make optional output argument)
      p_etsl=etsl;
      % also, assign etsl and the numerical code for behaviors to optional output
      % arguments
      varargout{1}=p_etsl;
      varargout{2}=AP.segmentType(:,1);

      % triggers with too short IEI - set them to 'too_brief' 
      % note: converting the length of an interval necessitates intv=0
      tmp2=cont2discrete(etsl(:,etslc.durCol),WP.osi*1e-3,'intv',0);
      etsl(tmp2<AP.ppSeg,etslc.tagCol)=briefV;
      logstr{end+1}=[int2str(length(find(tmp2<AP.ppSeg))) ' segments too brief & marked as such on behavioral event list'];
      disp(logstr{end});
      
      clear tmp* tthresh netsl petsl h fh_*
      % job following this one MUST be 'gen_bseq' 
      AP.job{1}='gen_bseq';
    
    % based on tsl from artifacts and behavioral data, generate list of segments to extract  
    % and show plot if requested
    case 'gen_bseq'       
      logstr{end+1}='---------- generating final list of usable segments..';
      disp(logstr{end});
      % save current etsl in temporary var (i for intermediate) for display
      i_etsl=etsl;
      % interweave etsls
      if ~isempty(af_etsl)
        etsl=etslindent(etsl,af_etsl);
        % keep a copy of etsl at this point. note that p_etsl is
        % overwritten here; further above it had been assigned to the etsl
        % devoid of artifact ts
        p_etsl=etsl;
        % alas, the whole business about too brief intervals has to be
        % repeated. See comment above.
        tmp2=cont2discrete(etsl(:,etslc.durCol),WP.osi*1e-3,'intv',0);
        tmp2=tmp2<AP.ppSeg & etsl(:,etslc.tagCol)~=artifV;
        etsl(tmp2,etslc.tagCol)=briefV;
        logstr{end+1}=[int2str(length(find(tmp2))) ' segments too brief & marked as such on behavioral event list'];
        disp(logstr{end});
      end
      % put out final etsl, if so desired
      varargout{3}=etsl;
      % etslindent adjusts duration column, so no need to recalculate here
      figure(WP.sumFigH); subplot('position',[.02 .85 .88 .125]);
      xlim=[0 WP.eor]/6e4;
      ylim=[.5 4.5];
      % time unit for plotting: min
      % attach any type of legal trigger to etsls used in current plot so last but one 
      % episode will be plotted, too
      r_etsl(end+1,[1 etslc.tagCol])=[WP.eoe badV];
      i_etsl(end+1,[1 etslc.tagCol])=[WP.eoe badV];        
      etsl(end+1,[1 etslc.tagCol])=[WP.eoe badV];
      axis([xlim ylim]);
      ph=repmat(nan,size(r_etsl,1)+size(i_etsl,1)+size(etsl,1)-1,1);
      % plot raw etsl
      for i=1:size(r_etsl,1)-1
        col=AP.segmentType{r_etsl(i,etslc.tagCol),3};
        ph(i)=patch([r_etsl(i:i+1,1); r_etsl([i+1,i],1)]/6e4,[3;3;2;2]+1.2,col);
        set(ph(i),'edgecolor','none');
      end
      text(xlim(2)*.5,4.2,'raw');                
      % intermediate 
      for j=1:size(i_etsl,1)-1
        col=AP.segmentType{i_etsl(j,etslc.tagCol),3};
        ph(i+j)=patch([i_etsl(j:j+1,1); i_etsl([j+1,j],1)]/6e4,[3;3;2;2],col);
        set(ph(i+j),'edgecolor','none');
      end
      text(xlim(2)*.5,3.1,'revised #1 (purged short segments & bad triggers)');
      % final
      for k=1:size(etsl,1)-1
        col=AP.segmentType{etsl(k,etslc.tagCol),3};
        ph(i+j+k)=patch([etsl(k:k+1,1); etsl([k+1,k],1)]/6e4,[3;3;2;2]-1.2,col);
        set(ph(i+j+k),'edgecolor','none');
      end
      text(xlim(2)*.5,1.9,'revised #2 (purged artifacts) ');
      niceyax(30);
      xlabel('time (min)');
      set(gca,'ytick',[],'box','on','xminortick','on');
      th=title(WP.figTitle);
      set(th,'fontsize',12,'fontweight','bold');
      % delete last triggers
      r_etsl(end,:)=[];      
      i_etsl(end,:)=[];      
      etsl(end,:)=[];
      % combo etsl is only needed for legend
      combo_etsl=[r_etsl; i_etsl; etsl];
      % commmand below does not work with fake last triggers in etsl
      [uet,tmp1]=unique(combo_etsl(:,etslc.tagCol));
      lh=legend(ph(tmp1),AP.segmentType(combo_etsl(tmp1,etslc.tagCol),1),4);
      set(lh,'position',[.91 .85 .07 .125]);
      % --- pie chart
      ep=unique(etsl(:,etslc.tagCol));
      col=cell(1);
      for i=1:length(ep)
        eptt(i)=sum(etsl(etsl(:,etslc.tagCol)==ep(i),etslc.durCol));
        col(i)=AP.segmentType(ep(i),3);          
      end
      subplot('position',[WP.xmarg .8/3*2+WP.ymarg .33/2-2*WP.xmarg .8/3-2*WP.ymarg]);
      % the handles 'pie' spits out are pointing to patches AND text alternately..
      ph=pie(eptt);
      for i=1:length(ep)
        set(ph(2*i-1),'facecolor',col{i},'edgecolor','none');
      end
      drawnow
      
      clear i j k r_etsl i_etsl ph uet tmp* combo* col ep* chld*
      AP.job(1)=[];

    case {'ovRawPlot'}
      % 'overview raw plot'
      logstr{end+1}='---------- generating overview plot of raw data..';
      disp(logstr{end});
      % fine line widths 
      labelscale('fontSz',8,'scaleFac',1.0,'lineW',.125,'markSz',4);
      fhr=mkfig(partJob); 
      rmouse_ovrawplot(p_etsl,partJob)
      delete(fhr)
      % reset graphics to what they were
      labelscale('fontSz',7,'scaleFac',1.0,'lineW',1.0,'markSz',4); 
      AP.job(1)=[];
      
    case {'rawPlot','csdPlot','chechao'}
      % all of the jobs below perform operations on short excerpts of raw
      % data, hence the preliminary work is combined
      switch partJob
        case 'rawPlot'
          logstr{end+1}='---------- plotting excerpts of raw data..';
          disp(logstr{end});
          % for plots of raw data we need fine line widths and not so fine
          % fonts
          labelscale('fontSz',10,'scaleFac',1.0,'lineW',.25,'markSz',4); 
          excLen=[30 5 .5];
        case 'csdPlot'
          logstr{end+1}='---------- plotting current source density of theta streams..'; 
          disp(logstr{end});
          % excerpt length should not exceed 10 s, otherwise not much will be discernible 
          % on plot
          excLen=[8 2 .5];
        case 'chechao' 
          logstr{end+1}='---------- plausibility check of channel order..';
          disp(logstr{end});
          % for a mosaic of tiny plots we need medium line widths and tiny fonts
          labelscale('fontSz',5,'scaleFac',1.0,'lineW',1.0,'markSz',8); 
          cmap=colormap(bone);
          % only first value is needed (keep others for compatibility reasons)
          excLen=[30 5 .5];
      end
      % take first decent stretch of data during exploratory behavior
      % without artifacts
      ix=find(etsl(:,etslc.tagCol)==explV & etsl(:,etslc.durCol)>=max(excLen*1000));
      if isempty(ix)
        parttitl='behavior undef';
        % pick any combination of segments devoid of bad signals/artifact
        tmpetsl=etsl;
        tmpetsl(etsl(:,etslc.tagCol)~=badV & etsl(:,etslc.tagCol)~=artifV,etslc.tagCol)=explV;
        tmp=~(diff(tmpetsl(:,etslc.tagCol)));
        tmpetsl([false; tmp],:)=[];
        tmp=diff([tmpetsl(:,etslc.tsCol); WP.eoe]);      
        tmpetsl(:,etslc.durCol)=tmp;
        ix=find(tmpetsl(:,etslc.tagCol)==explV & tmpetsl(:,etslc.durCol)>=max(excLen*1000));
        if isempty(ix)
          warning('job(s) rawPlot|csdPlot|chechao: no single stretch of data without artifacts or bad episodes could be found - using arbitrary one now');
          ix=find(tmpetsl(:,etslc.tagCol)==explV);
        end
        start=tmpetsl(ix(1),etslc.tsCol)*.001;
      else
        parttitl='exploring';
        start=etsl(ix(1),etslc.tsCol)*.001;        
      end
      % pick middle parts of segments
      start=start*[1 1 1]+cumsum([0 (excLen(1)-excLen(2))/2 (excLen(2)-excLen(3))/2]);
      stop=start+excLen;
      % cell array of strings for automatic adjustment of raw plot scale
      tmpscl=cell(1,AP.nCh);
      [tmpscl{:}]=deal(rawCh.dType);
      mkfig(partJob); 
      
      % job-specific parts:
      switch partJob
        case 'rawPlot'
          orient tall
          % load & plot all 2b analyzed
          for i=1:3
            subplot(3,1,i)
            if exist([dpath '\' abfFn '.mat'],'file')
              [d,si]=matDload([dpath '\' abfFn '.mat'],'start',start(i),'stop',stop(i),'channels',AP.rawChAnNm);        
            elseif exist([dpath '\' abfFn '.raw'],'file')
              d=rawload([dpath '\' abfFn '.raw'],AP.rawChAnNm,[start(i) stop(i)]*1000,rawFInfo);
              % put into array and convert to mV
              d=cat(2,d{:})/1000;
              tmp=rawfi2genfi(rawFInfo,'channels',AP.rawChMonNm);
              si=tmp.si;
            else
              [d,si]=abfload([dpath '\' abfFn '.abf'],'start',start(i),'stop',stop(i),'channels',AP.rawChAnNm);        
            end
            if DS.rawSignalInverted
              d=-1*d;
            end
            if i==1
              % assess range of amplitudes across channels
              dy=3*std(d(1:10:end));
            end
            pllplot(d,'si',si,'yscale',tmpscl,'spacing','fixed','dy',dy);
            if i==1,  
              titl=[abfFn ', ' parttitl ', ' sprintf('t=[%4.1f %4.1f]', [start(i) stop(i)]) ' s'];
            else
              titl=[parttitl ', ' sprintf('t=[%4.1f %4.1f]', [start(i) stop(i)]) ' s'];    
            end
            title(titl);      
          end
        case 'csdPlot'
          colormap(jet);
          % issue at least a warning if distances between electrodes are not uniform
          if length(unique(diff(cat(1,DS.rawCh{AP.LFPIdx,3}))))>1,
            logstr{end+1}='inter-electrode distances are not uniform';
            disp(logstr{end});
          end
          orient landscape
          nLev=20;
          for i=1:3
            % preallocate with nans (for missing channels)
            intvPts=cont2discrete([start(i) stop(i)]*1e6,WP.osi,'intv',1);
            d=repmat(nan,diff(intvPts)+1,AP.nAllLFPCh-2);
            for ci=1:nLFPCh
              d(:,AP.LFPccInd(ci))=strmread([AP.strmDir '\' rawCh(AP.LFPInd(ci)).thetaFn],'intv',intvPts,'verbose',0);
            end
            d=cursd(d);
            subplot(3,1,i),
            [c,cph]=contourf(discrete2cont((intvPts(1):intvPts(2)),WP.osi*1e-6),WP.elx(2:end-1),d',nLev);
            clear c;
            set(cph(:),'linestyle','none'); 
            cph=gca;
            % don't forget to flip axis
            set(cph,'ytick',WP.elx(2:end-1),'ydir','reverse');
            axis tight
            % line to indiciate principal channel
            line(get(gca,'xlim'),WP.elx(AP.LFPpcInd2)*[1 1],'color','k','linestyle','--','linewidth',1);
            if i==1,  title([abfFn ', ' parttitl]); end
            xlabel('time (s)');
            ylabel('electrode depth (mm)');
            % scale
            colorbar;
          end
        case 'chechao'
          % issue a warning if distances between chosen channels are not uniform
          if length(unique(diff(cat(1,DS.rawCh{AP.LFPIdx,3}))))>1,
            logstr{end+1}='inter-electrode distances are not uniform';
            disp(logstr{end});
          end
          orient landscape
          if exist([dpath '\' abfFn '.mat'],'file')
            [d,si]=matDload([dpath '\' abfFn '.mat'],'start',start(1),'stop',stop(1),'channels',AP.rawChAnNm(AP.LFPInd));        
          elseif exist([dpath '\' abfFn '.raw'],'file')
            d=rawload([dpath '\' abfFn '.raw'],AP.rawChAnNm,[start(1) stop(1)]*1000,rawFInfo);
            % put into array and convert to mV
            d=cat(2,d{:})/1000;
            tmp=rawfi2genfi(rawFInfo,'channels',AP.rawChMonNm);
            si=tmp.si;
          else
            [d,si]=abfload([dpath '\' abfFn '.abf'],'start',start(1),'stop',stop(1),'channels',AP.rawChAnNm(AP.LFPInd));        
          end
          if DS.rawSignalInverted
            d=-1*d;
          end
          % some data have nasty offsets, and we're not interested in delta, 
          % so filter at 5 Hz; also get rid of freqs > 40 Hz
          d=bafi(d,si,[5 40],'rs',40);
          % for this task with hippocampal data 100 ms should suffice as CC length..
          lag=cont2discrete(100,WP.osi*.001,'intv',1);
          % ..and also 100 ms the interval within which to spot a cow, ahem, peak
          cow=cont2discrete(99,WP.osi*.001,'intv',1);
          % cell arr of CCs
          c=cell(AP.nLFPCh);
          c(:)={repmat(nan,2*lag+1,1)};
          % matrix of peak corr coeffs
          cm=repmat(nan,AP.nLFPCh,AP.nLFPCh);
          % since the graphs have to be plotted sequentially anyways, do the XC within loop, channel
          % by channel
          for i1=1:AP.nLFPCh
            for i2=1:AP.nLFPCh
              [c{i1,i2},lags]=xcorr(d(:,i1),d(:,i2),lag,'coeff');
              tmpr=evdeal(c{i1,i2}(lag+1-cow:lag+1+cow),'idx','minmaxpeak');
              cm(i1,i2)=tmpr.maxPeak;
              if isnan(cm(i1,i2)), 
                cm(i1,i2)=-1;
              end
              subplot(AP.nLFPCh,AP.nLFPCh,(i1-1)*AP.nLFPCh+i2),
              % inflate plots
              rexy('xfac',1.2,'yfac',1.2);
              plot(lags,c{i1,i2},'r');
              set(gca,'ylim',[-.6 1.3]);
              grid on
              % set bg color of plot according to peak cc (mapping [-1 1] to [1 length of colormap])
              set(gca,'color',cmap(round(cm(i1,i2)*31)+32,:));
              title([AP.rawChAnNm{AP.LFPInd(i1)} ' vs ' AP.rawChAnNm{AP.LFPInd(i2)}]);              
            end
          end
      end
      drawnow  
      if ~isempty(AP.printas{1}), 
        for i=1:length(AP.printas)
          pa=AP.printas{i};
          if strfind(pa,'ps'), ext='.ps';
          elseif strfind(pa,'jpeg'), ext='.jpg';
          else ext='';
          end
          if strcmp(partJob,'chechao')
            % make background of plots visible
            set(gcf,'inverthard','off')   
          end
          print(pa,[WP.figName '_' partJob ext]);   
        end
      end
      % reset graphics to what they were
      labelscale('fontSz',7,'scaleFac',1.0,'lineW',1.0,'markSz',4); 
      set(gcf,'inverthard','on')   
      
      clear tmp* excLen ix partt* start stop d i titl pa ext nLev c cph i1 i2 cm cow lag lags
      AP.job(1)=[];
      
    case 'specgram'
      logstr{end+1}='---------- computing spectrogram..';
      disp(logstr{end});
      rmouse_specgram
      AP.job(1)=[];
      
    case 'filter&hilbert'
      logstr{end+1}='---------- filtering & Hilbert transform..';
      disp(logstr{end});
      rmouse_genstreams;
      AP.job(1)=[];             
      
    case 'gen_avDiff'      
      logstr{end+1}='---------- computing average raw slope..';
      disp(logstr{end});
      rmouse_genDiffTrace
      AP.job(1)=[];             
      
    case 'thetaPeaks'
      logstr{end+1}='---------- detecting peaks in theta streams..';
      disp(logstr{end});
      strmType='theta';
      rmouse_detstreampeak(rawCh,strmType);
      AP.job(1)=[];             
      
    case 'gammaEnvPeaks'
      logstr{end+1}='---------- detecting peaks in gammaEnv streams..';
      disp(logstr{end});
      strmType='gammaEnv';
      rmouse_detstreampeak(rawCh,strmType);
      AP.job(1)=[];             

    case 'segPrep'
      logstr{end+1}='---------- preparing segment-dependent computations..';
      disp(logstr{end});
      % remove all artifact, too brief etc. events from etsl
      etsl(etsl(:,etslc.tagCol)==badV,:)=[];
      etsl(etsl(:,etslc.tagCol)==briefV,:)=[];
      etsl(etsl(:,etslc.tagCol)==artifV,:)=[];
      etsl(etsl(:,etslc.tagCol)==undefV,:)=[];      
      % *** also, do without grooming because this behavior has so far only
      % taken up resources but never been interesting
      etsl(etsl(:,etslc.tagCol)==grooV,:)=[];            
      % create stop column
      etsl(:,stopCol)=etsl(:,etslc.tsCol)+etsl(:,etslc.durCol);
      % ** etsl in ticks **
      etsl_pts=etsl;
      etsl_pts(:,[etslc.tsCol stopCol])=cont2discrete(etsl(:,[etslc.tsCol stopCol]),WP.osi*.001,'intv',0);
      % duration in ticks (useful for preallocation)
      etsl_pts(:,etslc.durCol)=diff(etsl_pts(:,[etslc.tsCol stopCol]),1,2)+1;
      % row and col indices for CC comp - LFP channels only
      WP.nccix=[];
      tmpc=repmat((1:AP.nLFPCh),AP.nLFPCh,1);
      tmpr=rot90(tmpc,-1);
      for i=1:AP.nLFPCh
        WP.nccix=[WP.nccix; [diag(tmpr,i-1) diag(tmpc,i-1)]];
        WP.diagNccix{i}=[diag(tmpr,i-1) diag(tmpc,i-1)];
      end
      tmpr=[]; 
      % generate templates of 1. major CC results variables and 2. derived parameters
      % variables (such as mean, std etc.): both are cell arrays, the row and column indices 
      % coding for LFP channel indices in DS; each nonempty, non-nan cell holds the mean CC 
      % curve or derived parameter; empty cells are those which would contain redundant data
      % (CC(1,2) is the same as CC(2,1) inverted); cells with a single scalar nan correspond 
      % to channels not requested for analysis. Correspondingly for auto- and cross spectral 
      % densities. Since all of this preparatory stuff is not computationally intensive don't 
      % bother to check whether really all templates are needed
      % - cell content of raw CC
      tmpccc=repmat(nan,[2*AP.ccLagPts+1 1]);
      % - cell content of power spectra: the length of individual data segments is AP.ppSeg, 
      % which should be a power of 2. Thus, there is no need to zero-pad data segments and 
      % the size of the one-sided psd is known: AP.ppSeg/2:
      tmppsc=repmat(nan,[AP.ppSeg/2 1]);      
      WP.ccTemplate=cell(AP.nAllLFPCh);
      WP.ccTemplate(AP.trixie)={nan};
      WP.psTemplate=WP.ccTemplate;
      % no need to create psDerTemplate; WP.ccDerTemplate can be used
      WP.ccDerTemplate=WP.ccTemplate;
      for chInd1=1:AP.nLFPCh
        for chInd2=chInd1:AP.nLFPCh
          WP.ccTemplate{AP.LFPccInd(chInd1),AP.LFPccInd(chInd2)}=tmpccc;
          WP.psTemplate{AP.LFPccInd(chInd1),AP.LFPccInd(chInd2)}=tmppsc;
        end
      end
      if ~exist([AP.resPath '\' AP.resFn '.mat'],'file') && strcmpi(AP.saveMode,'append')
        disp('''append'' saving mode not applicable because results file not existent');
        AP.saveMode='replace';
      end
      switch AP.saveMode
        case 'replace'
          % generate results struct array - as many elements as AP.segmentType
          % has rows, and in exactly the same order
          % - some of the fields below are generated only because somewhere further down
          % their contents are checked (so they must exist) 
          for i=1:nSegmentType
            r(i).segmentType=AP.segmentType{i,1};
            r(i).ePts=[];
            r(i).iPts=[];        
            r(i).nePts=0;
            r(i).ni=0;
            r(i).rawPMn={[]};
            r(i).rawPStd={[]};        
            r(i).rawCohMn={[]};
            r(i).rawCohStd={[]};
            r(i).deCCMn={[]};
            r(i).deCCStd={[]};            
            r(i).thCCMn={[]};
            r(i).thCCStd={[]};            
            r(i).thLoeCCMn={[]};
            r(i).thLoeCCStd={[]};            
            r(i).thHieCCMn={[]};
            r(i).thHieCCStd={[]};            
            r(i).gaCCMn={[]};
            r(i).gaCCStd={[]};            
            r(i).gaNaCCMn={[]};
            r(i).gaNaCCStd={[]};            
            r(i).gaeCCMn={[]};
            r(i).gaeCCStd={[]};            
            r(i).thgaeCCMn=[];        
            r(i).thgaeCCStd=[];                    
            r(i).thgaeCCEnvMn=[];        
            r(i).thgaeCCEnvStd=[];                    
            r(i).detheCCMn=[];        
            r(i).detheCCStd=[];                    
            % segment start and stop pts
            idx=etsl_pts(:,etslc.tagCol)==i;
            if any(idx)
              % segment start and stop pts
              r(i).ePts=etsl_pts(idx,[etslc.tsCol stopCol]);
              % # of segments (of variable size)
              r(i).ne=size(r(i).ePts,1);
              % total # of pts IN EPISODES
              r(i).nePts=sum(etsl_pts(idx,etslc.durCol));
              % computational interval start and stop pts
              r(i).iPts=[];
              % ** note: conversion of time done here (by mkintrvls) is slightly different
              % from the way done further up by cont2discrete (in job 'gen_btsl', look for line
              % tmp2=cont2discrete(etsl(:,etslc.durCol),WP.osi*1e-3,'intv',0);)
              % For mkintervals to produce the intended intervals the right boundary of r(i).ePts 
              % must be incremented by one
              for j=1:size(r(i).ePts,1)
                tmp=mkintrvls(r(i).ePts(j,:)+[0 1],'resol',1,'ilen',AP.ppSeg,'olap',AP.dftOlapPts,'border','skip');
                r(i).iPts=cat(1,r(i).iPts,tmp);
              end
              % total # of COMPUTATIONAL INTERVALS 
              r(i).ni=size(r(i).iPts,1);
              % last and first points (in file) belonging to this segType
              r(i).lastPt=r(i).iPts(end,2);
              r(i).firstPt=r(i).iPts(1,1);
              % memory required for holding all channels in workspace at once
              % 8 bytes per data point expressed in MB  
              r(i).dmem=round(r(i).ePts(end,2)*AP.nCh*8/2^20);
            end
          end
        case 'append'
          % load previous results
          logstr{end+1}='---------- loading results from previous run..';
          disp(logstr{end});
          load([AP.resPath '\' AP.resFn]);          
          % to make sure that crucial parameters are the same between the
          % current and the previous run check field ePts of r
          % problem code: 0=all OK; 1=disparity explained below; 2=any other disparity (leading to program break)
          problem=zeros(nSegmentType,1);
          tmpr.ePts=[];
          for i=1:nSegmentType
            idx=etsl_pts(:,etslc.tagCol)==i;
            if any(idx)
              % segment start and stop pts
              tmpr(i).ePts=etsl_pts(idx,[etslc.tsCol stopCol]);
              problem(i)=2*double(~isequal(tmpr(i).ePts,r(i).ePts));
%               % -----------------------------------------------------------
%               % in the course of a bug fix in Sep 05 the computation of the
%               % stop column of etsl changed such that the intervals had to be
%               % shortened by exactly one point. However, at least for the beta3 KO
%               % project all *.mat results files had already been computed.
%               % Hence running any job with the option AP.saveMd='append' on
%               % one of these data sets caused problems here. The following
%               % lines correct for this, and only this, writing results variable r with the
%               % corrected .ePts fields 
%               if problem(i)
%                 if isequal(size(tmpr(i).ePts),size(r(i).ePts))
%                   if isequal(tmpr(i).ePts-r(i).ePts,repmat([0 1],size(r(i).ePts,1),1))
%                     problem(i)=1;
%                     % *** re-compute numbers
%                     % segment start and stop pts
%                     r(i).ePts=etsl_pts(idx,[etslc.tsCol stopCol]);
%                     % # of segments (of variable size)
%                     r(i).ne=size(r(i).ePts,1);
%                     % total # of pts IN EPISODES
%                     r(i).nePts=sum(etsl_pts(idx,etslc.durCol));
%                     % computational interval start and stop pts
%                     r(i).iPts=[];
%                     for j=1:size(r(i).ePts,1)
%                       tmp=mkintrvls(r(i).ePts(j,:)+[0 1],'resol',1,'ilen',AP.ppSeg,'olap',AP.dftOlapPts,'border','skip');
%                       r(i).iPts=cat(1,r(i).iPts,tmp);
%                     end
%                     % total # of COMPUTATIONAL INTERVALS
%                     r(i).ni=size(r(i).iPts,1);
%                     % last point (in file) belonging to this segType
%                     r(i).lastPt=r(i).iPts(end,2);
%                     % memory required for holding all channels in workspace at once
%                     % 8 bytes per data point expressed in MB
%                     r(i).dmem=round(r(i).ePts(end,2)*AP.nCh*8/2^20);
%                   end
%                 end
%               end
            end % if:any(idx)
          end % for:nSegmentType
          % now differentiate between situations
          if any(problem)
            if any(problem==2)
              logstr{end+1}=['data file: ' AP.resPath '\' DS.abfFn '; scoring file: ' AP.resPath '\' AP.bScoreFn ' - ''append'' mode for saving results is illegal because at least one key analysis parameter in AP was changed. Stopped analysis. You have to re-run it in ''replace'' mode.'];
              disp(logstr{end});
              warndlg(logstr{end});
              break
            else
              logstr{end+1}=['data file: ' AP.resPath '\' DS.abfFn '; scoring file: ' AP.resPath '\' AP.bScoreFn ' - owing to a recent bug fix (Sep 05) r.ePts and other timing information needs to be updated (=written to the main results file)'];
              disp(logstr{end});
              warndlg(logstr{end});
              % ** do it
              save([AP.resPath '\' AP.resFn],'r');
            end
          end
        otherwise
          error('illegal string in AP.saveMode');
      end
      clear idx tmp* problem
      AP.job(1)=[];                   
      
    case 'seg_rawSpecPow'
      logstr{end+1}='---------- computing spectral power of raw data..';
      disp(logstr{end});
      rmouse_cspecp(rawCh);
      r=orderfields(r);
      save([AP.resPath '\' AP.resFn],'r');
      AP.job(1)=[];                   
      
    case {'seg_gammaEnvSpecPow','seg_gammaNarrowEnvSpecPow'}
      if strfind(partJob,'Narrow')
        logstr{end+1}='---------- computing spectral power of narrow gamma envelope..';
        strmType='gammaNarrowEnv';
        STshort='gaNae';
      else
        logstr{end+1}='---------- computing spectral power of gamma envelope..';
        strmType='gammaEnv';
        STshort='gae';
      end
      disp(logstr{end});
      rmouse_cspecp_envelope(rawCh,strmType,STshort);
      r=orderfields(r);
      save([AP.resPath '\' AP.resFn],'r');
      AP.job(1)=[];                   
      
    case 'seg_deltaCC'
      logstr{end+1}='---------- computing delta CC..';
      disp(logstr{end});
      strmType='delta';
      % shortcut used as part of name of some variables
      STshort='de';
      rmouse_cc(rawCh,strmType,STshort);
      r=orderfields(r);
      save([AP.resPath '\' AP.resFn],'r');
      AP.job(1)=[];

    case 'seg_thetaCC'
      logstr{end+1}='---------- computing theta CC..';
      disp(logstr{end});
      strmType='theta';
      % shortcut used as part of name of some variables
      STshort='th';
      rmouse_cc(rawCh,strmType,STshort);
      r=orderfields(r);
      save([AP.resPath '\' AP.resFn],'r');
      AP.job(1)=[];
      
    case 'seg_thetaLoEnvCC'
      logstr{end+1}='---------- computing low-freq theta envelope CC..';
      disp(logstr{end});
      strmType='thetaLoEnv';
      % shortcut used as part of name of some variables
      STshort='thLoe';
      rmouse_cc(rawCh,strmType,STshort);
      r=orderfields(r);
      save([AP.resPath '\' AP.resFn],'r');
      AP.job(1)=[];
      
    case 'seg_thetaHiEnvCC'
      logstr{end+1}='---------- computing hi-freq theta envelope CC..';
      disp(logstr{end});
      strmType='thetaHiEnv';
      % shortcut used as part of name of some variables
      STshort='thHie';
      rmouse_cc(rawCh,strmType,STshort);
      r=orderfields(r);
      save([AP.resPath '\' AP.resFn],'r');
      AP.job(1)=[];

    case 'seg_gammaCC'
      logstr{end+1}='---------- computing gamma CC..';
      disp(logstr{end});
      strmType='gamma';
      % shortcut used as part of name of some variables
      STshort='ga';
      rmouse_cc(rawCh,strmType,STshort);
      r=orderfields(r);
      save([AP.resPath '\' AP.resFn],'r');
      AP.job(1)=[];
      
    case 'seg_gammaNarrowCC'
      logstr{end+1}='---------- computing narrow gamma CC..';
      disp(logstr{end});
      strmType='gammaNarrow';
      % shortcut used as part of name of some variables
      STshort='gaNa';
      rmouse_cc(rawCh,strmType,STshort);
      r=orderfields(r);
      save([AP.resPath '\' AP.resFn],'r');
      AP.job(1)=[];

    case 'seg_gammaEnvCC'
      logstr{end+1}='---------- computing gamma envelope CC..';
      disp(logstr{end});
      strmType='gammaEnv';
      % shortcut used as part of name of some variables
      STshort='gae';
      rmouse_cc(rawCh,strmType,STshort);
      r=orderfields(r);
      save([AP.resPath '\' AP.resFn],'r');
      AP.job(1)=[];
      
    case 'seg_gammaNarrowEnvCC'
      logstr{end+1}='---------- computing narrow gamma envelope CC..';
      disp(logstr{end});
      strmType='gammaNarrowEnv';
      % shortcut used as part of name of some variables
      STshort='gaNae';
      rmouse_cc(rawCh,strmType,STshort);
      r=orderfields(r);
      save([AP.resPath '\' AP.resFn],'r');
      AP.job(1)=[];

    case 'seg_thetaGammaEnvCC'
      if AP.ccNShuffle          
        logstr{end+1}='---------- computing CC theta - gamma envelope with shuffling..';
      else
        logstr{end+1}='---------- computing CC theta - gamma envelope..';
      end
      disp(logstr{end});
      ccType='thgae';
      rmouse_cc_intrasite(rawCh,ccType);
      r=orderfields(r);
      save([AP.resPath '\' AP.resFn],'r');
      AP.job(1)=[];
      
    case 'seg_deltaThetaEnvCC'
      if AP.ccNShuffle          
        logstr{end+1}='---------- computing CC delta - theta envelope with shuffling..';
      else
        logstr{end+1}='---------- computing CC delta - theta envelope..';
      end
      disp(logstr{end});
      ccType='dethe';
      rmouse_cc_intrasite(rawCh,ccType);
      r=orderfields(r);
      save([AP.resPath '\' AP.resFn],'r');
      AP.job(1)=[];

    case 'seg_rawGammaEnvCoh'
      logstr{end+1}='---------- computing coherence raw - gamma envelope ..';
      disp(logstr{end});
      cohType='rawgae';
      rmouse_coh_strms(cohType,rawCh)
      r=orderfields(r);
      save([AP.resPath '\' AP.resFn],'r');
      AP.job(1)=[];

    case 'seg_cc_ms2rad'
      logstr{end+1}='---------- converting theta-related phase lags to radians..';
      disp(logstr{end});
      rmouse_cc_ms2rad;
      
      rmouse_rcleanup;
      
      r=orderfields(r);
      save([AP.resPath '\' AP.resFn],'r');
      AP.job(1)=[];
      
    case 'tcPrincComp'
      logstr{end+1}='---------- computing & plotting principal components of select variables..';
      disp(logstr{end});
      rmouse_princcomp_tc;
      AP.job(1)=[];  
      
    case 'tcParCorr'
      logstr{end+1}='---------- computing & plotting correlations between select variables..';
      disp(logstr{end});
      rmouse_parcorr_tc;
      AP.job(1)=[];  

    case 'tcFig'
      logstr{end+1}='---------- producing time course plots of select variables..';
      disp(logstr{end});
      rmouse_p_tcsel;
      AP.job(1)=[];  

    case 'seg_thetaPeakReg'
      logstr{end+1}='---------- assessing regularity of theta based on variance of peak amplitude and inter-peak interval';
      disp(logstr{end});
      % stream type
      strmType='theta';
      STshort='th';
      rmouse_segoscillregularity(strmType,rawCh,STshort)
      r=orderfields(r);
      save([AP.resPath '\' AP.resFn],'r');
      AP.job(1)=[];  
      
    case 'seg_gammaEnvPeakReg'
      logstr{end+1}='---------- assessing regularity of gamma envelope based on variance of peak amplitude and inter-peak interval';
      disp(logstr{end});
      % stream type
      strmType='gammaEnv';
      STshort='gae';
      rmouse_segoscillregularity(strmType,rawCh,STshort)
      r=orderfields(r);
      save([AP.resPath '\' AP.resFn],'r');
      AP.job(1)=[];  

    case 'tcThetaPeak'
      logstr{end+1}='---------- producing blob (time course) plots of theta peaks..';
      disp(logstr{end});
      rmouse_p_tcthpeaks(p_etsl)
      AP.job(1)=[];  
      
    case 'sumFig'
      logstr{end+1}='---------- plotting summary of results..';
      disp(logstr{end});
      rmouse_p;
      AP.job(1)=[];  
      
    case {'callJob','seg_callJob'}
      eval([callJob ';']);
      if exist('r','var') && isstruct(r),
        r=orderfields(r);
        save([AP.resPath '\' AP.resFn],'r');
      end
      AP.job(1)=[];
      
    otherwise
      logstr{end+1}=['data file: ' AP.resPath '\' DS.abfFn '; scoring file: ' AP.resPath '\' AP.bScoreFn ' ---------- illegal subjob in AP.job - deleting from to-do list'];      
      warndlg(logstr{end});
      disp(logstr{end});
      AP.job(1)=[];      
  end
  done=isempty(AP.job);
end

% last act: close logfile
logstr={['********** analysis end: ' datestr(now)]};
disp(logstr{end});
diary off

% global variables must NOT BE DELETED! 
WP=[]; DS=[]; AP=[];




% ------------- trash ------------- trash ------------- trash ------------- trash 
% to convert from min.sec to ms: b=floor(a)*60000+(a-floor(a))*1e5
% to convert from ms to min.sec: a=floor(b/60000)+(b/60000-floor(b/60000))*.6;
% % display string and append it to text file
% function logb(fn,s)
% disp(s);
% fid=fopen(fn,'at');
% fprintf(fid,'%s \n',s);
% fclose(fid);
