function etsl=etslindent(los,win,varargin)
% ** function etsl=etslindent(los,win,varargin)
% makes one etsl 'indent' the other. Imagine two separate etsl describing
% different aspects of one and the same experimental period, e.g. one etsl
% lists behavioral events and the other neuronal events. The task is to
% combine the two etsl into one such that the events in one (the winner
% etsl) 'override' those in the other (the loser). That is, if a 'w' event
% intersects with a 'l' event, the 'l' event will be cropped, split up or
% even kicked out if the overlap is complete. The most graphical
% description is that of an indentation of one etsl into the other, hence
% the function name.
% ****  see etslconst.m for information on variables of type 'etsl' **** 
%
%                    >>> INPUT VARIABLES >>>
%
% NAME        TYPE/DEFAULT        DESCRIPTION
% los         etsl                the 'loser' etsl 
% win         etsl                the 'winner' etsl
%
%                    <<< OUTPUT VARIABLES <<<
%
% NAME        TYPE/DEFAULT    DESCRIPTION
% etls        etsl            result

verbose=0;
pvpmod(varargin);

etslconst;
if verbose, disp(['** ' mfilename]); end

wsz=size(win,1);
lsz=size(los,1);

% the column holding event START times
startcol=etslc.tsCol;
% temporary new column: event STOP times
stopcol=etslc.tagCol+1;
durcol=etslc.durCol;
tagcol=etslc.tagCol;
% checks
if any(win(1:end-1,startcol)+win(1:end-1,durcol)>=win(2:end,startcol))
  error('winner etsl contains overlapping or strictly adjacent events');
end
if any(los(1:end-1,startcol)+los(1:end-1,durcol)>los(2:end,startcol))
  error('loser etsl contains overlapping events');
end
% 'working' etsl: combination of both 
[etsl,idx]=sortrows([win; los],startcol);
clear win los;
% ** idx is the index to the 'winner' events within the new etsl: the loop
% below picks each of them, determines whether there is any overlap with
% other events (which can only be of 'loser' type) and adjusts these
% accordingly
[nix,idx]=intersect(idx,1:wsz);
% ** since the loop progresses along winner events, in all loop runs
% following the first one there is no need to check loser events which had
% already been found not to overlap with the previous winner event. i1ind
% points to the first loser event that can possibly have an overlap.
i1ind=1;
% create temporary new column  (start + duration)
etsl(:,stopcol)=etsl(:,startcol)+etsl(:,durcol);

% one of the 4 possible 'overlaps' between winner and loser (case 4, which
% cuts a loser event in 2 pieces) necessitates inserting events in the
% etsl. Since it would be very time-consuming to enlarge the size of etsl
% every time that happens preallocate as many lines as there are winner
% events. Setting their value to nan prevents them being detected as valid
% events
etsl(wsz+lsz+1:wsz+lsz+1+wsz,:)=nan;
nane=repmat(nan,size(etsl(1,:)));

for i=1:length(idx)
  % [i1 & i4] point to events with any overlap with the current win event:
  % the leftmost: loser events terminating later than current winner event
  % begins
  i1=etsl(i1ind:end,stopcol)>=etsl(idx(i),startcol);
  % the rightmost: loser events starting before current winner event
  % terminates
  i4=etsl(i1ind:end,startcol)<=etsl(idx(i),stopcol);
  % reference for cidx: whole etsl
  cidx=find(i1 & i4)+i1ind-1;
  % i1 & i4 will always contain the current winner event, so let's kick it
  % out  
  cidx=setdiff(cidx,idx(i));
  if ~isempty(cidx)
    i2=etsl(cidx,stopcol)<=etsl(idx(i),stopcol);
    i3=etsl(cidx,startcol)>=etsl(idx(i),startcol);
    % case 1:
    %       .......
    %  ......     .......  loser
    %    
    %           .......
    %   .........     ...  winner
    ii= i2 & ~i3;
    if any(ii)
      if length(find(ii))>1, error('internal:case1 should be max 1 event'); end
      % stop time of loser event = start time of winner event
      etsl(cidx(ii),stopcol)=etsl(idx(i),startcol);    
    end

    % case 2:
    %           .......
    %   .........     ...  loser
    %    
    %       .......
    %  ......     .......  winner
    ii= ~i2 & i3;
    if any(ii)
      if length(find(ii))>1, error('internal:case2 should be max 1 event'); end
      % start time of loser event = stop time of winner event
      etsl(cidx(ii),startcol)=etsl(idx(i),stopcol);    
    end

    % case 3:
    %         ...
    %   .......  .......  loser
    %    
    %      .........
    %  .....       .....  winner
    ii= i2 & i3;
    if any(ii)
      % delete all losers
      etsl(cidx(ii),:)=nan;    
    end

    % case 4: winner splits loser in 2
    % *** this is the only case generating new events; since the whole code
    % relies on the etsl being properly sorted, it MUST come last and
    % insert the new event in the right place 
    %       .........
    %  ......       .....  loser
    %    
    %          ...
    %  ......... ........  winner
    ii= ~i2 & ~i3;
    if any(ii)
      if length(find(ii))>1, error('internal:case4 should be max 1 event'); end
      % create new loser event
      e=nane;
      e([startcol stopcol tagcol])=[etsl(idx(i),stopcol)  etsl(cidx(ii),stopcol)  etsl(cidx(ii),tagcol)];
      % adjust stop time of original loser event
      etsl(cidx(ii),stopcol)=etsl(idx(i),startcol);
      % intercalate new loser event (must come right after current winner
      % event)
      etsl(idx(i)+2:end,:)=etsl(idx(i)+1:end-1,:);
      etsl(idx(i)+1,:)=e;
      % don't forget to adjust idx..
      idx(i:end)=idx(i:end)+1;
    end
  end
  % re-calculate i1 and adjust i1ind for the next loop run
  i1=etsl(i1ind:end,stopcol)>etsl(idx(i),startcol);
  if any(i1), i1ind=i1ind+min(find(i1))-1; end
end

% finale:
% 1. delete nan rows 
etsl(isnan(etsl(:,startcol)),:)=[];
% 2. check whether start and stop times are OK (=no overlaps of events)
% ** importante: detection of overlap must be eps-tolerant
tmp=makecol(rot90(etsl(:,[stopcol startcol])));
if any(diff(tmp)<tmp(1:end-1)*(-eps)), error('internal: overlapping events'); end
% 3. compute duration
etsl(:,durcol)=diff(etsl(:,[startcol stopcol]),1,2);
% 4. delete temporary column
etsl(:,stopcol)=[];
