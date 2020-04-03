function varargout=combine_bv(varargin)
% ** function varargout=combine_bv(varargin)
% - combines individual rmouse data sets and computes and plots relative amount 
%   of time spent in different behaviors (the time course of it)
% - needs ANPAR and DSET (=concatenation of individual and matching AP and DS)
%   as global variable
% - the (optional) input argument mt is an array of points in time corresponding 
%   to each experiment's synchronization point. These could be e.g. the points when 
%   some drug was administered. The unit is MINUTES from the start of recording of 
%   the first file. If not specified, everything will be aligned according to start points
%   of files in first row of DSET/ANPAR
%
%                         >>> INPUT VARIABLES >>>
%
% NAME          TYPE/DEFAULT          DESCRIPTION
% mt            1d array              see above; mt must have as many elements as
%                                     ANPAR/DSET have columns (experiments)
% resol         scalar, .5            time resolution in seconds for the plot
%                                     showing time evolution of behavior. Too
%                                     high values will result in incorrect
%                                     depiction of short segments, too low
%                                     values are computationally inefficient


global ANPAR DSET AP DS

mt=[];
resol=.5;
pvpmod(varargin)

[n1 n2 n3]=size(ANPAR);
[dn1 dn2 dn3]=size(DSET);

% if mt not specified, align everything to start points of files in first row
if isempty(mt)
  % @ check this
  mt=zeros(1,n2);
else
  % convert mt to seconds for internal use
  mt=mt*60;
end
mtn=length(mt);

if ~isequal([n1 n2 n3],[dn1 dn2 dn3]), error('DSET and ANPAR must have equal dimension'); end
if ~numel(DSET), error('DSET and ANPAR are empty or, more likely, not declared global in the base workspace'); end
if mtn~=n2, error('the number of elements in mt must be equal to the number of columns in DSET and ANPAR'); end
 
% choose behaviors to be plotted (legal value of AP.segmentType)
% ** last one must be 'undef'
behav={'grooming','immobile','exploring','undef'};
nBehav=length(behav);
dispFormatString=repmat('%8.3f',1,nBehav);

rmouse_ini;
etslconst;
stgp;

% generate behavior-specific index into segTypeGlobP
for bi=1:nBehav
  stgpix(bi)=strmatch(behav{bi},segTypeGlobP(:,1),'exact');
end
if length(stgpix)~=nBehav, error('check behavior/stgp'); end

close all;
labelscale('fontSz',10,'scaleFac',1.0,'lineW',1.0,'markSz',6); 


% -------- PART I: collection of data

% the compound etsl
cetsl=[];
% string containing file name and proportion of time spent in behaviors 
bPropStr=cell(n1*n2,1);
bPropStr(:)={['empty_slot: ' repmat(' nan',1,nBehav)]};
% loop over experiments
fct=0;
globalBTime=zeros(1,nBehav);
for ci=1:n2
  % 'local' (for each experiment) etls
  etsl=[];
  % loop over concentrations
  for ri=1:n1
    AP=ANPAR(ri,ci);
    DS=DSET(ri,ci);
    fct=fct+1;
    % since experiments may be composed of a variable number of files, some elements of
    % ANPAR and DSET may be empty structures (test one vital field of each to ascertain)
    if ~isempty(DS.abfFn) & ~isempty(AP.resFn)
      % if paths do not contain a drive letter, pre-pend WP.rootPath
      if isempty(strfind(DS.dpath,':')), DS.dpath=[WP.rootPath DS.dpath]; end
      % extract si and timing information from abf file -
      % if matfile exists, pick it instead of abf file
      if exist([DS.dpath '\' DS.abfFn '.mat'],'file')
        load([DS.dpath '\' DS.abfFn '.mat'],'abfi');
      elseif exist([DS.dpath '\' DS.abfFn '.abf'],'file')
        [nix,nix2,abfi]=abfload([DS.dpath '\' DS.abfFn '.abf'],'info');
        % abfi=abfinfo([DS.dpath '\' DS.abfFn '.abf'],'verbose',0);
        % put out start and stop times in seconds from midnight
        disp(['recording start||end: ' sprintf('%6.3f || %6.3f',abfi.recTime) ' s from midnight']);
      else
        error([DS.dpath '\' DS.abfFn ' does not exist'])
      end
      % offset in seconds
      tOffs=abfi.lFileStartTime;
      if ri==1
        % this is the offset to be subtracted from all time stamps of the current
        % experiment
        syncTOffs=tOffs+mt(ci);
      end
      % run rmouse with output argument
      AP.job={'gen_btsl'};
      % the trouble is that rmouse.m clears global variables internally, so we need a
      % local copy to keep
      % tmpAP=AP;
      % @ consider using a default AP
      [e,rmSegType]=rmouse;
      % AP=tmpAP;
      % convert etsl to seconds and add offset
      e(:,[etslc.tsCol, etslc.durCol])=e(:,[etslc.tsCol, etslc.durCol])*.001;
      e(:,etslc.tsCol)=e(:,etslc.tsCol)+tOffs;
      % indices to behaviors
      bix={[]};
      for bi=1:nBehav
        bix{bi}=e(:,etslc.tagCol)==strmatch(behav{bi},rmSegType(:,1),'exact');
      end
      % combine all behaviors different from the ones requested and tag them as 'undef'
      % (=include them in the last category)
      bix{end}=bix{end} | ~any(cat(2,bix{:}),2);
      % e(~any(cat(2,bix{:}),2),etslc.tagCol)=nan;
      % total duration of recording (s)
      recDur=sum(e(:,2));
      % in this loop
      % (i) change numerical codes of behavior to those of segTypeGlobalP
      % (ii) compute proportion of time spent in certain behaviors 
      bProp=[];
      for bi=1:nBehav
        e(bix{bi},etslc.tagCol)=stgpix(bi);
        globalBTime(bi)=globalBTime(bi)+sum(e(bix{bi},etslc.durCol));
        bProp(bi)=sum(e(bix{bi},etslc.durCol))/recDur;
      end
      % collect relative amounts of time spent in each behavior
      bPropStr{fct}=[sprintf('%30s', [DSET(ri,ci).dpath '\' DSET(ri,ci).abfFn ': ']) num2str(bProp,dispFormatString) ];
      % concatenate etsls
      etsl=cat(1,etsl,e);
      
    end % if isstruct(AP) & isstruct(DS)
  end % for ri=1:n1
  % this is essentially the same test as ~isempty(DS.abfFn) & ~isempty(AP.resFn)
  if ~isempty(etsl)
    % align all events of current experiment to sync point
    etsl(:,etslc.tsCol)=etsl(:,etslc.tsCol)-syncTOffs;
    % concatenate
    cetsl=cat(1,cetsl,etsl);
  end
end

% dump info on screen
disp('*** behaviors analyzed:');
disp(behav);
disp('*** normalized amounts of time spent in these behaviors:')
for ii=1:n1*n2
  disp(bPropStr{ii});
end

% *** if one or more of the behaviors in variable behav do not occur 
% in the file(s) analyzed, the output of etsliefreq will not represent 
% these behaviors (=ief and nief will contain less than nBehav columns). 
% As a remedy, insert columns of zeroes into ief and nief for the 
% missing behaviors. This is somewhat stupid but provides an output which
% facilitates combining results from separate runs of combine_bv, should
% the need ever arise. So, 
% 1. generate inst ev freq
[t,ief,nief]=etsliefreq(cetsl,resol);
% 2. the columns of ief and nief are sorted according to the numerical 
% value of the tag, so stgpix has to be sorted
[stgpix,ix]=sort(stgpix);
% 3. reorder globalBTime accordingly and use it as an index into behaviors
% present in files analyzed
globalBTime=globalBTime(ix);
% 4. finally, rearrange/insert columns
iefTemplate=zeros(length(t),nBehav);
iefTemplate(:,find(globalBTime))=ief;
ief=iefTemplate;
iefTemplate(:,find(globalBTime))=nief;
nief=iefTemplate;

% spit out time series, if requested
varargout{1}=t;
varargout{2}=ief;
varargout{3}=nief;
varargout{4}=segTypeGlobP(stgpix,1);

clear e etsl cetsl bix iefTemplate

% ---------- PART II: plot -------------------------------------------------
ftag='behavTimeCourse';
figha=findobj('tag',ftag);
if isempty(figha), figha=figure;
else  figure(figha);
end
tmpScrSz=get(0,'Screensize')*.8+.1*rand;
tmpScrSz([1 2])=tmpScrSz([1 2])+tmpScrSz([3 4])*.10;
set(figha,'position',tmpScrSz,'tag',ftag,'name',ftag,...
  'color',[.9 .9 1],'numbertitle','off');
clf;
orient landscape;
labelscale('fontSz',10,'scaleFac',1,'lineW',1,'markSz',6); 

% time axis in minutes
subplot(2,1,1)
ah=area(t/60,ief);
ct=0;
for bi=find(globalBTime)
  ct=ct+1;
  set(ah(ct),'facecolor',segTypeGlobP{stgpix(bi),2})
end
subplot(2,1,2)
ah=area(t/60,nief);
ct=0;
for bi=find(globalBTime)
  ct=ct+1;
  set(ah(ct),'facecolor',segTypeGlobP{stgpix(bi),2})
end
xlabel('time (min)');

% % marker
% if rvi==1 & ~isempty(mt)
%   mph=plot(mt,repmat(co(1)+.9*diff(co),size(mt)),'rv');
%   set(mph,'markersize',16,'markerfacecolor','y');
% end
% % expand plots a little
% rexy('xfac',1,'yfac',1.5);

