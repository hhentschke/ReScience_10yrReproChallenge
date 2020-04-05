function [cutout,isCutout,ntslc]=tsl2exc(d,si,tslc,varargin)
% ** function [cutout,isCutout,ntslc]=tsl2exc(d,si,tslc,varargin)
%    extracts fixed-length cutouts from continuous raw data d.
%    Cutouts are listed in time stamp list tsl.
%
%                    >>> INPUT VARIABLES >>>
%
% NAME           TYPE/DEFAULT              DESCRIPTION
% d              1d or 2d arr              sampled data, time runs along columns
% si             a) scalar                 sampling interval in us
%                b) string 'idx'           means that ALL timing information 
%                                            is specified in points (as opposed to ms)
%                                            and will also be put out as such
% tslc           CELL ARRAY of 1d arr      time stamp list (tsl) 
% .........................................................................
%       NOTE:
%         - if only one tsl is given, this tsl is used to produce cutouts 
%           from all columns of d
%         - if as many tsl are given as columns in d, cutouts will be 
%           produced from matching pairs tsl/column in d
%         - any other combination of numbers of tsl/columns of d is illegal
% .........................................................................
% win            2-element arr, [-1 1]     cutout window relative to time stamp (ms or pts)
% deadT          scalar, 0                 dead time (ms or pts)
%
%                    <<< OUTPUT VARIABLES <<<
%
% NAME           TYPE/DEFAULT               DESCRIPTION
% cutout          3d-array                   cutouts: <time> <cutout#> <channel#>
%                   OR
%                cell array of 2d arr       cutouts (along columns), from 1 channel per element
%                                             if cutout is not completely contained in d,
%                                             it will not be put out
% isCutout       cell array of 1d arr       an index into the time stamp lists in tslc;
%                                             1=cutout could be produced for time stamp
%                                             0=cutout could NOT be produced (was on the border of 
%                                             or outside of sweep
% ntslc          cell array of 1d arr       time stamp lists devoid of events not respecting 
%                                             borders/dead time


% improvements:

win=[1 1];
deadT=0;
verbose=0;
pvpmod(varargin);

if numel(win)~=2
  error('input par ''win'' must have exactly two elements')
end

if verbose, disp(['**** ' mfilename ':']); end
nTsl=length(tslc);

if ischar(si)
  if ~strcmpi(si,'idx')
    error('input variable si must be either a scalar or string ''idx''');
  end
  winPts=win;
  deadTPts=deadT;
else
  % ** convert to ticks (easier to calculate with) **
  winPts=cont2discrete(win,si*.001,'intv',1);
  deadTPts=cont2discrete(deadT,si*.001,'intv',1);
  for g=1:nTsl
    tslc{g}=cont2discrete(tslc{g},si*.001);
  end
end
preTrigPts=winPts(1);
postTrigPts=winPts(2);

coNPts=diff([preTrigPts postTrigPts])+1;
if coNPts<1, error('check pre/posttrigger'); end
if deadTPts<0, error('check dead time'); end

[dn1 dn2]=size(d);
if nTsl==1
  % preallocate cutouts
  cutout=zeros(coNPts,nTsl,dn2);
  isMultipleTsl=false;
elseif nTsl>1 & nTsl==dn2
  % initialize cutout
  cutout=cell(1,dn2);
  isMultipleTsl=true;
else  
  error([mfilename ': the dimensions of d and tslc must be either [<any> 1] or [n n]']);
end;

v1=[];
isCutout=cell(1,dn2);
ntslc=cell(1,dn2);
for g=1:nTsl
  % time stamps expressed in points (take only first column in case an etls
  % was provided inadvertently)
  tsPts=tslc{g}(:,1);
  % number of time stamps
  nTs=length(tsPts);
  % temporary variable containing cutouts
  v1=[];
  % cutouts - one per column
  if ~isempty(tsPts)
    % first locate events fully contained in sweep by inspecting border points 
    v1=repmat([preTrigPts postTrigPts],nTs,1)+repmat(tsPts,1,2);
    isCutout{g}=((v1(:,1)>=1) & (v1(:,2)<=dn1));
    % delete others from list
    tsPts(~isCutout{g})=[];
    % dead time:
    if deadTPts>0
      tsPts=tsldeadt(tsPts,deadTPts);
    end
    %  NEW number of time stamps..
    nTs=length(tsPts);
    %  ..and tsl
    ntslc{g}=tsPts;
    if nTs
      v1=repmat([preTrigPts:postTrigPts]',1,nTs)+repmat(tsPts',coNPts,1);
      if isMultipleTsl
        v1=d(v1,g);
        % <time> <cutout#> 
        v1=reshape(v1,coNPts,nTs);
      else
        v1=d(v1,:);
        % <time> <cutout#> <channel#>
        v1=reshape(v1,coNPts,nTs,dn2);
      end
    else
      v1=zeros(coNPts,0,nTsl);
    end
  end
  if iscell(cutout)
    cutout{g}=v1;
  else
    % this can only be the case if only one tsl was given and d is an array
    % so the loop will be executed only once
    cutout=v1;
  end
end

if ~ischar(si)
  % no need to check whether si==idx
  for g=1:nTsl
    if ~isempty(ntslc{g})
      ntslc{g}=discrete2cont(ntslc{g},si*.001);
    end
  end
end

