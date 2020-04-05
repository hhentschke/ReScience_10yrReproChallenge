function rcat(behav,rv,varargin)
% ** function rcat(behav,rv,varargin)
% - concatenates selected fields of results variable r distributed over files
% - needs ANPAR and DSET as global vars: these must be 1-column struct arrays, the 
%   elements consisting of individual and matching AP and DS 
% - expects (empty) global variable R 
%
%                         >>> INPUT VARIABLES >>>
% NAME           TYPE/DEFAULT          DESCRIPTION
% behav          cell array of chars   behaviors that should make it in the final output array
%                                      any (combination) of {'immobile','exploring'}
% rv             cell array of chars   any (combination) of fields of r containing
%                                      segment-wise computed parameters, like 'thCCPeak'
%                     
%                         <<< OUTPUT VARIABLES <<<
%
% NAME             TYPE/DEFAULT     DESCRIPTION
% R (global)       struct           contains fields as specified in input var rv,
%                                   all of them 3D arrays: 
%                                             n by n by t
%                                                  OR
%                                             1 by n by t,
%                                   n being the number of channels, time along slices,
%                                   + an extra field TB: 
%                                   1st column is time in seconds from start of recording
%                                   of first file (midpoint of each interval)
%                                   2nd colum is behavioral tag

global ANPAR DSET R

pvpmod(varargin);

[n1 n2 n3]=size(ANPAR);
[dn1 dn2 dn3]=size(DSET);
if ~isequal([n1 n2 n3],[dn1 dn2 dn3]), error('DSET and ANPAR must have equal dimension'); end
if n2>1 | n3>1, error('DSET and ANPAR must have no more than 1 column'); end

nrv=length(rv);
% make sure that major output var R is empty
R=[];
% R: contains fields as specified in rv
for rvi=1:nrv
  eval(['R.' rv{rvi} '=[];']);
end
% variable collecting global time & behavType
TB=[];

rmouse_ini;
% this generates variable segTypeGlobP containing global settings for plots of behaviors
stgp;
bcode=[];
for i=1:length(behav)
  bcode=[bcode strmatch(behav{i},segTypeGlobP(:,1),'exact')];
end
if length(bcode)~=length(behav)
  error('unlikely, but true: input par ''behav'' has entries not found in segTypeGlobP');
end

% -------- PART I: set up time & master templates and compute slice indices for each
% segment within the 3D master vars
% loop over data sets for one experiment
for ii=1:length(ANPAR)
  AP=ANPAR(ii);
  DS=DSET(ii);
  % if paths do not contain a drive letter, pre-pend WP.rootPath
  if isempty(strfind(DS.dpath,':')), DS.dpath=[WP.rootPath DS.dpath]; end
  if isempty(strfind(AP.resPath,':')), AP.resPath=[WP.rootPath AP.resPath]; end    
  % extract si and timing information from abf file - if matfile exists, pick it instead of abf file
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
  % original sampling interval
  osi=abfi.si;
  % load results var..
  load([AP.resPath '\' AP.resFn],'r');
  % channel business
  rawCh=rmouse_chan;
  % index to last row of current TB
  lix=size(TB,1);
  % time offset of current file to add to current time (seconds)
  tOffs=abfi.lFileStartTime;
  if ii==1
    % time offset of first file (seconds) which should be the reference time t=0
    masterTOffs=tOffs;
  end
  %   if lix>0
  %     tOffs=abfi.lFileStartTime+discrete2cont(round(AP.ppSeg/2),osi*.001)/1e3;
  %   else 
  %     tOffs=abfi.lFileStartTime;
  %   end

  % loop over behaviors, collect parameters
  tb=[];
  for ri=1:length(r)
    lbi=strmatch(r(ri).segmentType,behav,'exact');
    if ~isempty(r(ri).iPts) & ~isempty(lbi)
      % time axis (seconds) for current file, reference being 0:00 a.m. 
      % of the day the first file in DSET/ANPAR was recorded
      t=discrete2cont(round(mean(r(ri).iPts,2)),osi*.001)/1e3+tOffs;
      nt=length(t);
      % 'global' behavioral code
      if isfinite(DS.conc) & DS.conc~= 0 & DS.conc~= -1 
        lbi=strmatch([behav{lbi} '_drug'],segTypeGlobP(:,1),'exact');
      else
        lbi=strmatch(behav{lbi},segTypeGlobP(:,1),'exact');
      end
      % local (=for current file) time & behav code
      tb=cat(1,tb,[t repmat(lbi,nt,1)]);
    end
  end
  
  % for each of the behaviors, tbix knows the slice index of each segment's results within the
  % master array for the current file
  [tb,tbix]=sortrows(tb,1);
  tb(:,3)=(1:size(tb,1))';
  tbix=tb(tbix,3);
  tb(:,3)=[];
  % ** this line makes the assumption that files were put in chronological order **
  TB=[TB; tb];

  % -------- PART II: collection of data
  for rvi=1:nrv
    Y=[];
    for ri=1:length(r)
      lbi=strmatch(r(ri).segmentType(:,1),behav);
      if ~isempty(r(ri).iPts) & ~isempty(lbi)
        % assign field contents to generic variable y, put in proper shape
        eval(['y=r(ri).' rv{rvi} ';']);
        % select extraction method depending on results var chosen (in essence,
        % differentiate between cell arrays and arrays)
        switch rv{rvi}
          case {'thgaeCCPeak','thgaeCCPeakT'}
            if ~isempty(y)
              % easy: permute
              y=permute(y,[3 2 1]);
              Y=cat(3,Y,y);
            else
              y=[];
            end
          otherwise
            if length(diag(y))>1
              % transfer values from cell array into array - has to be done elementwise 
              % because there are empty cells & we're dealing with 3D arrays
              yy=repmat(nan,[nAllLFPCh nAllLFPCh r(ri).ni]);              
              for rw=1:nAllLFPCh
                for cl=rw:nAllLFPCh
                  yy(rw,cl,:)=permute(y{rw,cl},[1 3 2]);
                end
              end
              Y=cat(3,Y,yy);
              clear yy
            else
              y=[];
            end
        end % switch
        if isempty(y)
          error(['data to be extracted do not exist']);
        end
      end % if:~isempty r.ipts & behav requested
    end % for:behaviors
    % put in right order
    Y=Y(:,:,tbix);
    % collect values
    eval(['R.' rv{rvi} '=cat(3,R.' rv{rvi} ',Y);']);
    y=[];
  end % for:rv
end

% subtract time offset
TB(:,1)=TB(:,1)-masterTOffs;
% make field of R
R.tb=TB;
% by default, attach AP and DS of first file
R.AP=ANPAR(1);
R.DS=DSET(1);