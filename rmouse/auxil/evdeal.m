function [r,d]=evdeal(d,si,job,varargin)
% ** function [r,d]=evdeal(d,si,job,varargin)
% extracts commonly useful parameters from multi-episode or multi-channel 
% time series data, like maxima, minima, peak amplitudes, etc.
%
%                         >>> INPUT VARIABLES >>>
% NAME             TYPE/DEFAULT             DESCRIPTION
% d                2d or 3d array           the waveforms (along columns) 
% si               a) scalar                sampling interval in us
%                  b) string 'idx'          means that all timing information in
%                                            results structure will be given in
%                                            points (as opposed to ms)
%                                            ** NOTE: you may get unexpected results
%                                            when calculating slopes!
% pretreat         cell array of strings,   determines how to pre-cook the data
%                    'halign'                align horizontally according to peak (includes zero-padding;
%                                            may triple memory requirement) - NOT YET IMPLEMENTED -
%                    'valign'                align vertically (=base line subtraction;
%                                            base line = leading 1/20 th of points)
% job              cell array of strings,   determines parameters to extract (if empty, nothing will be done)
%                    'all'                   all of the following
%                    'minmax'                min, max & max-min
%                    'minmaxslope'           min, max & max-min of SLOPE
%                    'minmaxpeak'            min and max peak (if no peak found, nan is returned)
%                    'allpeaks'              all pos and neg peaks, sorted (if no peak found, nan is returned)
%                    'closepeak'             pos and neg peak closest to point in time tP (see below) 
%                    'area'                  integral of abs(curve)
%                    'fwidth'                width at fraction 'rH' of max height of peak, see below 
%                                            **NOTE : fractional width will be determined from POSITIVE,  
%                                            POSITIVE-GOING peak of choice (see below), assuming baseline=0. 
%                    'pow'                   total power (spectral estimate)
%                    'rms'                   root mean square
%                    'pc'                    principal components (up to 8) - NOT YET IMPLEMENTED -
% plateauOK        scalar, 1                qualifier for any of the '...peak' jobs: if zero,
%                                            plateaus will not be recognized as peaks
%                                            (realize that in sampled data with poor
%                                            resolution of amplitude peaks will more often
%                                            than not be misrepresented as plateaus)
% tP               scalar or array, 0       qualifier for job 'closepeak'; must be given in ms or points
%                                            depending on choice of si; if not a scalar, number of 
%                                            columns must correspond to number of columns in d
% pTy              string, 'maxPeak'        qualifier for job 'fwidth'; must be either 'maxPeak' or 'tlPosPeak'
%                                            this choice determines the type of peak
%                                            **NOTE: choice of job 'fwidth' implies that either of jobs 
%                                            'minmaxpeak' or 'closepeak' will be run - make sure that parameters  
%                                            for the latter are set correctly 
% rH               scalar or array, 0.5     qualifier for job 'fwidth'; must be in ]0 1[ 
% verbose          scalar, 0                if nonzero you will be informed about details 
%                                           of the calculations
%
%                         <<< OUTPUT VARIABLES <<<
%
% NAME             TYPE/DEFAULT            DESCRIPTION
% d               same as d                modified (according to pretreatment) data
%                                            if no pretreatment had been specified, an empty array will be given 
% r                structure               the results; individual fields are 1d or 2d
%                                          arrays with as many elements as waveforms in d; 
%                                          names correspond to the job name/should be self-explanatory
%                   
% ** NOTE: if si is a scalar ALL TIME VARIABLES ARE ms, and time is mapped 
% such that the first data point in d corresponds to t=0 in real time. 
% That is, if the maximum of a trace happens to fall on the first point,
% the resultant .maxPeakT will be 0 (and not si)
% 

% **** NOTE: there are a few dependencies of jobs:
% - most jobs need diffd
% - 'fwidth' needs r.maxPeakT as ticks, not as ms 
% ---> the way these dependencies are taken into account relies on the order of jobs,
% which in the current version is controlled by the position of the code. Should this
% be changed (by implementing a while loop which 'eats away' each partjob (as in e.g.
% rmouse.m) this logic has to be reevaluated

% ----- default values & varargin -----
verbose=0;
pretreat={};
tP=0;
rH=.5;
pTy='maxPeak';
plateauOK=1;
pvpmod(varargin);
if verbose, disp(['**** ' mfilename ':']); end;

% -------------- preliminaries ----------------------------------------------------
% this may happen if variables are accidentally given as strings
if ~iscell(pretreat), pretreat={pretreat}; end
if ~iscell(job), job={job}; end  
% if si is string..
if ischar(si)
  if strcmpi(si,'idx')
    sieqpts=1;
    tP_pts=tP;
  else
    error('input variable si must be either a scalar or string ''idx''');
  end
else
  sieqpts=0;
  tP_pts=cont2discrete(tP,si*.001,'intv',0);
end
if rH<=0 | rH>=1,  error('rH must be e ]0 1['); end

md=[];
r.min=[];
r.minT=[];
r.max=[];
r.maxT=[];
r.totampl=[];
r.slopemin=[];
r.slopeminT=[];
r.slopemax=[];
r.slopemaxT=[];
r.minPeak=[];
r.minPeakT=[];
r.maxPeak=[];
r.maxPeakT=[];
% tl for time-locked
r.tlPosPeak=[];
r.tlPosPeakT=[];
r.tlNegPeak=[];
r.tlNegPeakT=[];
r.fwidth=[];
% attention, here cum cell arrays
r.posPeak={[]};
r.posPeakT={[]};
r.negPeak={[]};
r.negPeakT={[]};
r.area=[];
r.totpow=[];
r.rms=[];
r.pc=[];

if ~isempty(strmatch('all',job,'exact'))
  job={'minmax','minmaxslope','minmaxpeak','allpeaks','closepeak','area','fwidth','pow','rms','pc'};
end  
% if half width is one of the jobs, include either of jobs needed (if not in list)
% order of jobs in 'job' does not matter, it is imposed below
ji1=isempty(strmatch('fwidth',job,'exact'));
ji2=isempty(strmatch('minmaxpeak',job,'exact'));
ji3=isempty(strmatch('closepeak',job,'exact'));
if ~ji1
  if strmatch('maxPeak',pTy) & ji2
    job{end+1}='minmaxpeak';
  elseif strmatch('tlPosPeak',pTy) & ji3
    job{end+1}='closepeak';    
  end
end
[n1 n2 n3]=size(d);
[pn1 pn2 pn3]=size(tP);
if pn1+pn2+pn3~=3 & ~isequal([n2 n3],[pn2 pn3])
  error('check input par tP');
end

% -------------- pretreatment ----------------------------------------------------
% vertical alignment MUST come first
sidx=strmatch('valign',pretreat,'exact');
if ~isempty(sidx)
  if verbose, disp('---> vertical alignment (base line subtraction)'); end
  % by default, subtract mean of leading 1/20th of data
  blPts=round(n1/20);
  d=d-repmat(mean(d(1:blPts,:,:),1),[n1 1 1]);
  pretreat(sidx)=[];
end

sidx=strmatch('halign',pretreat,'exact');
if ~isempty(sidx)
  if verbose, disp('to come soon: ---> horizontal alignment (peaks)'); end
  % for this, we need the maximum, so lets put it in r and thus possibly prevent redundant calculation further down
  [r.max tmpMaxT]=max(d);
  if sieqpts, r.maxT=tmpMaxT;
  else r.maxT=discrete2cont(tmpMaxT,si*.001,'intv',0); end

  % DONT CONFUSE POINTS AND MS FOR CALCULATION OF SHIFT
  
  pretreat(sidx)=[];
end

% if anything is left in pretreat, we have a problem
if ~isempty(pretreat)
  warning(['one of pretreatment methods not recognized & ignored: ' pretreat ' (check spelling)']);
end

% -------------- parameter calculations ----------------------------------------------------
sidx=strmatch('minmaxslope',job,'exact');
if ~isempty(sidx)
  if verbose, disp('---> determining min and max of SLOPES'); end
  diffd=diff(d);
  [tmpSMin tmpSMinT]=min(diffd);
  [tmpSMax tmpSMaxT]=max(diffd);
  if sieqpts
    % unit: uV/tick
    r.slopemin=tmpSMin;
    r.slopemax=tmpSMax;
    % unit: ms
    r.slopeminT=tmpSMinT;
    r.slopemaxT=tmpSMaxT;
  else
    % unit: uV/ms
    r.slopemin=1e3*tmpSMin/si;
    r.slopemax=1e3*tmpSMax/si;
    % unit: ms
    r.slopeminT=discrete2cont(tmpSMinT,si*.001,'intv',0);
    r.slopemaxT=discrete2cont(tmpSMaxT,si*.001,'intv',0);
  end
  job(sidx)=[];  
end

sidx=strmatch('minmax',job,'exact');
if ~isempty(sidx)
  if verbose, disp('---> determining min, max & max-min'); end
  % if r.max and maxT exist, assume they had been calculated OK above
  if isempty(r.max) | isempty(r.maxT)
    [r.max tmpMaxT]=max(d);
    if sieqpts, r.maxT=tmpMaxT;
    else r.maxT=discrete2cont(tmpMaxT,si*.001,'intv',0); 
    end
  end
  [r.min tmpMinT]=min(d);
  if sieqpts, r.minT=tmpMinT;
  else r.minT=discrete2cont(tmpMinT,si*.001,'intv',0); 
  end
  r.totampl=r.max-r.min;
  job(sidx)=[];
end

sidx=strmatch('minmaxpeak',job,'exact');
if ~isempty(sidx)
  if verbose, disp('---> determining min & max PEAKS'); end
  if ~exist('diffd','var'),
    diffd=diff(d);
  end
  % preallocation
  r.maxPeak=repmat(nan,[1 n2 n3]); 
  r.maxPeakT=r.maxPeak;
  r.minPeak=r.maxPeak;  
  r.minPeakT=r.maxPeak;  
  for j=1:n3
    for i=1:n2
      % tmp1: potential positive-going peaks = -1, potential negative-going peaks = 1
      tmp1=diff((diffd(:,i,j)>0),1,1);
      % positive-going peaks
      pgpi=find(tmp1==-1)+1;
      % tmp1 now: the reverse
      tmp1=diff((-1*diffd(:,i,j)>0),1,1);
      % negative-going peaks
      ngpi=find(tmp1==-1)+1;
      if ~plateauOK
        % accept only peaks in the strict sense: a point flanked by two points less in
        % absolute amplitude (otherwise it would be a point marking the beginning or
        % end of a plateau)
        if ~isempty(pgpi) 
          pgpi=pgpi(all(repmat(d(pgpi,i,j),1,2)>[d(pgpi-1,i,j) d(pgpi+1,i,j)],2));
        end
        if ~isempty(ngpi) 
          ngpi=ngpi(all(repmat(d(ngpi,i,j),1,2)<[d(ngpi-1,i,j) d(ngpi+1,i,j)],2));
        end
      end
      if ~isempty(pgpi)
        [r.maxPeak(1,i,j),ix]=max(d(pgpi,i,j));
        if sieqpts, r.maxPeakT(1,i,j)=pgpi(ix);
        else
          r.maxPeakT(1,i,j)=discrete2cont(pgpi(ix),si*.001,'intv',0);
        end
      end
      if ~isempty(ngpi)
        [r.minPeak(1,i,j),ix]=min(d(ngpi,i,j));
        if sieqpts, r.minPeakT(1,i,j)=ngpi(ix);
        else r.minPeakT(1,i,j)=discrete2cont(ngpi(ix),si*.001,'intv',0);
        end
      end
    end
  end
  job(sidx)=[];
end
  
sidx=strmatch('closepeak',job,'exact');
if ~isempty(sidx)
  if verbose, disp('---> determining max & min time-locked peaks'); end
  if ~exist('diffd','var'),
    diffd=diff(d);
  end
  % preallocation
  r.tlPosPeak=repmat(nan,[1 n2 n3]);  
  r.tlPosPeakT=r.tlPosPeak;
  r.tlNegPeak=r.tlPosPeak; 
  r.tlNegPeakT=r.tlPosPeak;
  for j=1:n3
    for i=1:n2
      % tmp1: potential positive-going peaks = -1, potential negative-going peaks = 1
      tmp1=diff((diffd(:,i,j)>0),1,1);
      % positive-going peaks
      pgpi=find(tmp1==-1)+1;
      % tmp1 now: the reverse
      tmp1=diff((-1*diffd(:,i,j)>0),1,1);
      % negative-going peaks
      ngpi=find(tmp1==-1)+1;
      if ~plateauOK
        % accept only peaks in the strict sense: a point flanked by two points less in
        % absolute amplitude (otherwise it would be a point marking the beginning or
        % end of a plateau)
        if ~isempty(pgpi) 
          pgpi=pgpi(all(repmat(d(pgpi,i,j),1,2)>[d(pgpi-1,i,j) d(pgpi+1,i,j)],2));
        end
        if ~isempty(ngpi) 
          ngpi=ngpi(all(repmat(d(ngpi,i,j),1,2)<[d(ngpi-1,i,j) d(ngpi+1,i,j)],2));
        end
      end
      if ~isempty(pgpi)
        [tmpm,ix]=min(abs(pgpi-tP_pts));
        r.tlPosPeak(1,i,j)=d(pgpi(ix),i,j);
        if sieqpts, r.tlPosPeakT(1,i,j)=pgpi(ix);
        else r.tlPosPeakT(1,i,j)=discrete2cont(pgpi(ix),si*.001,'intv',0);
        end
      end
      if ~isempty(ngpi)
        [tmpm,ix]=min(abs(ngpi-tP_pts));
        r.tlNegPeak(1,i,j)=d(ngpi(ix),i,j);
        if sieqpts, r.tlNegPeakT(1,i,j)=ngpi(ix);
        else r.tlNegPeakT(1,i,j)=discrete2cont(ngpi(ix),si*.001,'intv',0);
        end
      end
    end
  end
  clear ngpi pgpi tmp1 tmpm
  job(sidx)=[];    
end

sidx=strmatch('allpeaks',job,'exact');
if ~isempty(sidx)
  if verbose, disp('---> determining ALL PEAKS'); end
  if ~exist('diffd','var'),
    diffd=diff(d);
  end
  % preallocation
  r.posPeak=cell([1 n2 n3]); r.posPeak(:)={nan};
  r.posPeakT=r.posPeak;
  r.negPeak=r.posPeak;
  r.negPeakT=r.posPeak;
  for j=1:n3
    for i=1:n2
      % tmp1: potential positive-going peaks = -1, potential negative-going peaks = 1
      tmp1=diff((diffd(:,i,j)>0),1,1);
      % positive-going peaks
      pgpi=find(tmp1==-1)+1;
      % tmp1 now: the reverse
      tmp1=diff((-1*diffd(:,i,j)>0),1,1);
      % negative-going peaks
      ngpi=find(tmp1==-1)+1;
      if ~plateauOK
        % accept only peaks in the strict sense: a point flanked by two points less in
        % absolute amplitude (otherwise it would be a point marking the beginning or
        % end of a plateau)
        if ~isempty(pgpi) 
          pgpi=pgpi(all(repmat(d(pgpi,i,j),1,2)>[d(pgpi-1,i,j) d(pgpi+1,i,j)],2));
        end
        if ~isempty(ngpi) 
          ngpi=ngpi(all(repmat(d(ngpi,i,j),1,2)<[d(ngpi-1,i,j) d(ngpi+1,i,j)],2));
        end
      end
      if ~isempty(pgpi)
        r.posPeak{1,i,j}=d(pgpi,i,j);
        if sieqpts, r.posPeakT{1,i,j}=pgpi;
        else r.posPeakT{1,i,j}=discrete2cont(pgpi,si*.001,'intv',0);
        end
      end
      if ~isempty(ngpi)
        r.negPeak{1,i,j}=d(ngpi,i,j);
        if sieqpts, r.negPeakT{1,i,j}=ngpi;
        else r.negPeakT{1,i,j}=discrete2cont(ngpi,si*.001,'intv',0);
        end
      end
    end
  end
  clear ngpi pgpi tmp1
  job(sidx)=[];    
end

sidx=strmatch('fwidth',job,'exact');
if ~isempty(sidx)
  if verbose, disp('---> determining half widths of peaks'); end
  if ~exist('diffd','var'),
    diffd=diff(d);
  end
  % preparations:
  % assign peak amplitudes and times depending on choice of pTy
  eval(['pA=r.' pTy ';']);
  eval(['pT=r.' pTy 'T;']);  
  % it's a drag, but if times are given in real units, we have to reconvert..
  if ~sieqpts
    pT=cont2discrete(pT,si*.001,'intv',0);
  end
  % preallocation
  r.fwidth=nan*zeros([1 n2 n3]); 
  for j=1:n3
    for i=1:n2
      % only if max pos peak exists width at fractional height can be determined
      if isfinite(pA(1,i,j))
        % tmp1: potential left intersection = 1, potential right intersection = -1
        tmp1=diff((d(:,i,j)>pA(1,i,j)*rH),1,1);
        % left intersection = 1
        ngpi=find(tmp1==1)+1;
        % right intersection = -1
        pgpi=find(tmp1==-1)+1;
        % the intersections, if any, must flank the peak. If either intersection 
        % does not exist, forget about it
        if ~isempty(ngpi) & ~isempty(pgpi)
          ngpi=ngpi(ngpi<pT(1,i,j));
          pgpi=pgpi(pgpi>pT(1,i,j));
        end
        % since ngpi and pgpi are (still) sorted, we can just pick the right ones
        if ~isempty(ngpi) & ~isempty(pgpi)
          if sieqpts, r.fwidth(1,i,j)=pgpi(1)-ngpi(end);
          else r.fwidth(1,i,j)=discrete2cont(pgpi(1)-ngpi(end),si*.001,'intv',1);
          end
        end
      end
    end
  end
  clear ngpi pgpi tmp1 pA pT
  job(sidx)=[];    
end

sidx=strmatch('area',job,'exact');
if ~isempty(sidx)
  if verbose, disp('---> determining area'); end
  r.area=sum(abs(d));
  % r.area=sum((d));  
  job(sidx)=[];  
end

sidx=strmatch('pow',job,'exact');
if ~isempty(sidx)
  if verbose, disp('---> determining power'); end
  winPts=n1;
  win=hanning(winPts);
  % for length of window, find the next larger power of 2 
  nfft=2^(ceil(log2(winPts)));
  for pidx2=1:n2
    for pidx3=1:n3
      [P,F]=periodogram(d(:,pidx2,pidx3),win,nfft);
      r.totpow(1,pidx2,pidx3)=sum(P);  
    end
  end
  job(sidx)=[];  
end

sidx=strmatch('rms',job,'exact');
if ~isempty(sidx)
  if verbose, disp('---> computing root mean square'); end
  r.rms=sqrt(mean(d.^2));
  job(sidx)=[];  
end

sidx=strmatch('pc',job,'exact');
if ~isempty(sidx)
  if verbose, disp('to come soon: ---> calculation of principal components'); end
  
  job(sidx)=[];  
end

% if anything is left in job, we have yet another problem
if ~isempty(job)
  warning(['one of the jobs not recognized & ignored: ' job ' (check spelling)']);
end
