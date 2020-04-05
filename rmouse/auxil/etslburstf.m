function [etsl,atsl,silentEtsl,varargout]=etslburstf(tsl,maxIEI,varargin)
% ** function [etsl,atsl,silentEtsl,varargout]=etslburstf(tsl,maxIEI,varargin)
%     combines time stamps fulfilling a burst criterion and puts out the
%     burst time stamp list in two different formats (etsl & atsl). See
%     etslconst.m, primer_evltsl.rtf and definition_atsl.rtf for
%     information about variables of type 'tsl', 'etsl' and 'atsl'.
%     Short description of how bursts are found:
%     i) all groups of time stamps separated by less than maxIEI from
%     each other are initially considered bursts, and the first events
%     within these groups are the potential burst starts
%     ii) if maxIEI_init is a positive number SMALLER than maxIEI, events
%     within above-mentioned groups will be further inspected: the earliest
%     pair of events separated by less than or equal to maxIEI_init will be
%     defined as the burst start (more precisely, the first event within
%     this pair will be the burst start. In other words, we define the
%     first instance of upshooting event rates as burst starts). If none is
%     found, the group of events in question will not be considered a
%     burst.
%     iii) if minSilentPerDur is a positive number a burst will only be
%     accepted if maximally maxSilentPerNEv events occur in this pre-burst
%     interval
%     iv) if minNEvPerBurst is a positive number >= 2 only bursts consisting of
%     at least minNEvPerBurst events will be accepted
%     v) if minNEvPerBurst is 1, single events satisfying condition i)
%     (distance to preceding event >= maxIEI) but not condition ii) will
%     also be considered bursts
%
%                    >>> INPUT VARIABLES >>>
% NAME            TYPE/DEFAULT     DESCRIPTION
% tsl             single col arr   the time stamp list
% maxIEI          scalar           maximum inter event interval
% maxIEI_init     scalar, maxIEI   maximum INITIAL inter event interval -
%                                  all time stamps separated by less than
%                                  this value will be considered burst
%                                  starters. Set to Nan if this criterion
%                                  shall not apply (see above)
% minSilentPerDur scalar, NaN      minimal duration (ms) of pre-burst
%                                  silent period
% maxSilentPerNEv scalar, NaN      maximal number of events allowed in
%                                  pre-burst silent period
% minNEvPerBurst  scalar, NaN      minimal number of events a burst must
%                                  consist of - watch out, if set to 1
%                                  single events will be considered bursts!
% startTs         char, 'ev'      'ev': burst start is set to first event
%                                  in bursts 
%                                  'iei': burst start is set to center of
%                                  initial inter-event-interval of bursts
%                                  (except when bursts may consist of just
%                                  one event)
% recLen          scalar, NaN      length of recording (ms)
%
%                    <<< OUTPUT VARIABLES <<<
% NAME            TYPE/DEFAULT    DESCRIPTION
% etsl            extended time stamp list
% atsl            advanced time stamp list
% silentEtsl      extended time stamp list of silent periods
% varargout{1}    struct with several statistics on bursts and silent
%                  periods

% to do
% - safety net for recs with only one burst or silent per

% ---- output args:
% etslburstf will only produce two of the columns of an etsl, so put out
% zero-by-two empty etsls in case tsl is empty or no burst is detected (but
% see further below for special treatment of atsl)
etsl=zeros(0,2);
atsl=zeros(0,2);
silentEtsl=zeros(0,2);
% default values of stats
if nargout>3
  stats.recTime=nan;
  stats.fractionEvInBurst=nan;
  stats.relTimeInBurst=nan;
  stats.mnNEvPerBurst=nan;
  stats.stdNEvPerBurst=nan;
  stats.mnIntraBurstEvRate=nan;
  stats.stdIntraBurstEvRate=nan;
  stats.mnBurstLen=nan;
  stats.stdBurstLen=nan;
  stats.burstRate=nan;
  stats.mnSilentPerLen=nan;
  stats.stdSilentPerLen=nan;
  varargout{1}=stats;
end

% ---- set defaults
etslconst;
maxIEI_init=maxIEI;
minSilentPerDur=nan;
maxSilentPerNEv=NaN;
minNEvPerBurst=NaN;
startTs='ev';
recLen=NaN;
pvpmod(varargin);

% ------------------------- checks of input -----------------------------------
if isempty(tsl)
  warning('time stamp list is empty - returning empty results');
  return
end

if isfinite(maxIEI)
  if maxIEI<=0.0
    error('maxIEI must be strictly positive');
  end
else
  error('maxIEI must be strictly positive and finite');
end

if isfinite(maxIEI_init)
  if maxIEI_init<=0.0
    error('maxIEI_init must be strictly positive');
  end
  if maxIEI<maxIEI_init
    error('maxIEI must be greater than maxIEI_init');
  end
else
  error('maxIEI_init must be strictly positive and finite');
end

if sum(isnan([minSilentPerDur maxSilentPerNEv]))==1
  error('input parameters ''minSilentPerDur'' and ''maxSilentPerNEv'' must BOTH be NaN or non-NaN');
end
if maxSilentPerNEv<0
  error('maxSilentPerNEv must be >=0');
end
if minNEvPerBurst<1
  error('minNEvPerBurst must be >=1');
elseif minNEvPerBurst==1
  % this is a somewhat funny request, so issue a warning
  warning([mfilename ':singleEvBu'],'considering single events as bursts');
end
if ~ismember(startTs,{'ev','iei'})
  error('startTs must be either of ''ev'' or ''iei''')
end
  
% last check:
if isfinite(recLen) 
  if recLen<tsl(end)
    error('''recLen'' is less than the last time stamp');
  end
else
  warning('recording length not specified - setting to last event in tsl (which may bias statistics');
  recLen=tsl(end);
end

% in case an etsl is put in retain just the ts column
[n1, n2]=size(tsl);
if n2>1
  tsl=tsl(:,etslc.tsCol);
end

% ---------------- main job ----------------------------------------
% atsl has to be initialized here (third col is needed for spike count)
atsl=[tsl zeros(n1,2)];
% the time differences between adjacent time stamps
delta=diff(tsl);
nDelta=length(delta);
% index to time stamps separated by at least maxIEI from precursor, plus
% the first event in the list(in dubio pro reo...). In other words, these
% are the first events occurring after the burst as defined thus far.
% They are thus the potential burst starters.
tmpIdx=[1; find(delta>maxIEI)+1];
if isempty(tmpIdx)
  warning('not a single event can be tagged ''beginning of a burst'' - putting out empty vars');
  % important to kill the third col of atsl
  atsl(:,3)=[];
else
  if maxIEI_init<maxIEI
    % within the groups of spikes delimited by potential burst starters
    % pick the first interval <= maxIEI_init, i.e. the first in a series of
    % events in quick succession: these are the true burst starters
    buInitIdx=nan*tmpIdx;
    nBuInit=length(buInitIdx);
    tmpIdx2=find(delta<=maxIEI_init);
    for h=1:nBuInit
      if h<nBuInit
        % tmpIdx(h):tmpIdx(h+1)-1 is the index to all events embraced by
        % potential burst start events (in a [ [ manner). tmpIdx2 points
        % to events IEI < maxIEI_init ('true' burst starters, that is,
        % those events initiating a burst of more than one event). There
        % may be more than one of the latter in each group of spikes, so
        % we have to find events satisfying both criteria
        localBuInitIdx=intersect((tmpIdx(h):tmpIdx(h+1)-1),tmpIdx2);
      else
        localBuInitIdx=intersect((tmpIdx(h):nDelta),tmpIdx2);
      end
      if ~isempty(localBuInitIdx)
        % pick the first 'true' burst starter (this relies on the indices
        % being sorted by intersect)
        buInitIdx(h)=localBuInitIdx(1);
      end
    end
    % nans represent single events or groups of loosely spaced events
    % without a proper burst start as defined above, so kill them
    buInitIdx(isnan(buInitIdx))=[];
    if isempty(buInitIdx)
      etsl=zeros(0,2);
    end
  else
    buInitIdx=tmpIdx;
  end %if:maxIEI_init<maxIEI

  % if single events shall be considered bursts, merge with tmpIdx
  if minNEvPerBurst==1
    buInitIdx=union(buInitIdx,tmpIdx);
  end
  nBuInit=length(buInitIdx);

  % next criterion: maximal number of events in pre-burst interval
  if ~isnan(minSilentPerDur)
    % time interval to (maxSilentPerNEv+1)-th neighbor
    % ** note that this does NOT yet exclude events which are less than
    % minSilentPerDur ms away from beginning of recording! **
    deltaT_neighN=diff([[repmat(-inf,maxSilentPerNEv+1,1); tsl(1:end-maxSilentPerNEv-1)] tsl],1,2);
    goodIx=deltaT_neighN(buInitIdx)>=minSilentPerDur;
    buInitIdx=buInitIdx(goodIx);
    nBuInit=length(buInitIdx);
  end
  clear tmpI* good*
  % ----- find events belonging to bursts
  % initialize both tsls with three columns: the third column holds the
  % number of events per burst, which is possibly an exclusion criterion
  etsl=nan(nBuInit,3);
  for i=1:nBuInit
    if minNEvPerBurst>1 && strcmpi(startTs,'iei')
      % set beginning of burst to center of first inter-event-interval
      etsl(i,etslc.tsCol)=sum(tsl(buInitIdx(i)+[0 1]))/2;
    else
      % set beginning of burst to first spike in burst
      etsl(i,etslc.tsCol)=tsl(buInitIdx(i));
    end
    % running variable for loop
    g=0;
    % duration is defined as the time span between beginnings of
    % events, so if a burst consists of only one event its duration is
    % zero
    dur=0;
    while (buInitIdx(i)+g<=nDelta && delta(buInitIdx(i)+g)<=maxIEI)
      dur=dur+delta(buInitIdx(i)+g);
      g=g+1;
    end
    etsl(i,etslc.durCol)=dur;
    % let's misuse the etslc.amplCol for the number of events partaking
    % in each burst 
    etsl(i,etslc.amplCol)=g+1;
    % now that both duration and number of events/burst are known fill
    % the slots in atsl
    atsl(buInitIdx(i),[2 3])=[dur g];
  end % for i=1:nBuInit

  % ----- further exclusion criterion: number of events per burst (by
  % definition they cannot be less than 1 at this point, checking makes
  % sense only for minNEvPerBurst>1)
  if minNEvPerBurst>1
    badIx=etsl(:,etslc.amplCol)<minNEvPerBurst;
    etsl(badIx,:)=[];
    % excluding bursts in an atsl amounts to setting the duration (and
    % spike count) column to zero
    atsl(buInitIdx(badIx),[2 3])=0;
    % buInitIdx is still needed, so track all changes
    buInitIdx(badIx)=[];
    nBuInit=length(buInitIdx);
  end

  % catch emptied lists
  if isempty(etsl)
    warning('etslburstf:emptiedEtsl','no burst left in list after application of exclusion criteria - putting out empty lists');
  else

    % -------------------------------------------------------------------
    % Now take care of bursts and silent periods at the beginning and end
    % of the recording. They must be gotten rid of because they are -
    % except for extremely rare cases - incomplete (the recording started
    % somewhere between beginning and end) and including them would bias
    % the statistics.
    % Another important aspect to consider: bursts and silent periods
    % alternate by definition. Therefore, in order to get unbiased
    % statistics on their frequency and the relative time the system spends
    % in either state we have to regard the interval between ALIKE
    % TRANSITIONS as relevant. Either from the first transition sp->bu to
    % the last such transition or from the first bu->sp the last such
    % transition.
    % -------------------------------------------------------------------

    % Accomplish the above with the help of a silentEtsl:
    % - start of silent periods=end of bursts (first to last-1)
    silentEtsl(1:nBuInit-1,etslc.tsCol)=sum(etsl(1:end-1,[etslc.tsCol etslc.durCol]),2);
    % - durations: up to start of following burst (2nd to last)
    silentEtsl(:,etslc.durCol)=diff([silentEtsl(:,etslc.tsCol) etsl(2:end,etslc.tsCol)],1,2);

    % the silent periods etsl is 'clean' at this point in the sense that its
    % first and last entry are bordered by the end and beginning,
    % respectively, of confirmed bursts (and not by the beginning and end,
    % respectively, of the recording, which are more likely than not in the
    % middle of either event).
    % So, now check the bursts themselves.
    % BEGINNING:
    % - if the first burst is < minSilentPerDur away from the beginning kick
    % it (which is conservative for cases in which more than zero spikes in
    % the preceding silent period are allowed)
    if etsl(1,etslc.tsCol)<minSilentPerDur
      etsl(1,:)=[];
      atsl(buInitIdx(1),[2 3])=0;
      buInitIdx(1)=[];
      nBuInit=nBuInit-1;
    end
    % - if the true recording length is specified as an input arg test
    % whether the end of the last bu is less than maxIEI away from the end
    % of the recording. If that is so, kick it
    if isfinite(recLen) && sum(etsl(end,[etslc.tsCol etslc.durCol]))+maxIEI>recLen
      etsl(end,:)=[];
      atsl(buInitIdx(end),[2 3])=0;
      buInitIdx(end)=[];
      nBuInit=nBuInit-1;
    end
    % At this point, both etsls are 'clean' in the sense mentioned above.
    % The next thing to do, for statistics, is to compute the proper time
    % interval for rate measures (that is, the time interval between the
    % first and the last alike transition).

    % variable output arg: all kinds of stats
    if nargout>3
      % checking wich kind of event comes first and computing .recTime
      % makes sense only if there is at least one burst and one silent
      % period
      if size(etsl,1)>0 && size(silentEtsl,1)>0
        % if the first event not affected by border problems is a burst...
        if etsl(1,etslc.tsCol)<silentEtsl(1,etslc.tsCol)
          % ...time starts ticking with a transition
          %   silent period (unknown length) -> burst
          % and lasts to the end of the last silent period (no matter whether a
          % valid bursts follows or not)
          stats.recTime=sum(silentEtsl(end,[etslc.tsCol etslc.durCol]),2)-etsl(1,etslc.tsCol);
        else
          % ...otherwise, time starts ticking with a transition
          % burst (unknown length) -> silentPer
          stats.recTime=sum(etsl(end,[etslc.tsCol etslc.durCol]),2)-silentEtsl(1,etslc.tsCol);
        end
      end
      % - fraction of events partaking in bursts (etslc.amplCol contains
      % the number of events in each burst)
      stats.fractionEvInBurst=sum(etsl(:,etslc.amplCol),1)/n1;
      % - relative time spent in bursts
      stats.relTimeInBurst=sum(etsl(:,etslc.durCol),1)/stats.recTime;
      % - average number of spx/burst
      stats.mnNEvPerBurst=mean(etsl(:,etslc.amplCol),1);
      stats.stdNEvPerBurst=std(etsl(:,etslc.amplCol),0,1);
      % - average intra-burst firing rate (Hz)
      tmp=(etsl(:,etslc.amplCol)./etsl(:,etslc.durCol))*1000;
      stats.mnIntraBurstEvRate=mean(tmp,1);
      stats.stdIntraBurstEvRate=std(tmp,0,1);
      % - burst length (mean+-std)
      stats.mnBurstLen=mean(etsl(:,etslc.durCol));
      stats.stdBurstLen=std(etsl(:,etslc.durCol));
      % - burst rate (identical to rate of silent periods)
      stats.burstRate=min(size(silentEtsl,1),nBuInit)/(stats.recTime/1000);
      % - length of silent periods (mean+-std)
      stats.mnSilentPerLen=mean(silentEtsl(:,etslc.durCol),1);
      stats.stdSilentPerLen=std(silentEtsl(:,etslc.durCol),0,1);
      varargout{1}=stats;
    end
    % make sure we retain only .tsCol and .durCol
    etsl(:,etslc.amplCol)=[];
    atsl(:,3)=[];
    % ** convert time to seconds in atsl only **
    atsl=atsl/1000;
  end % if:isempty(etsl)
end % if:isempty(tmpIdx) (burst starters found)
