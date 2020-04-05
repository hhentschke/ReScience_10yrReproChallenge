function tsl=tcd(d,si,thresh,varargin)
% ** function tsl=tcd(d,si,thresh,varargin)
% 'threshold crossing detector' - a threshold-based spike detector
% producing time stamp lists
%
%                    >>> INPUT VARIABLES >>>
% NAME       TYPE/DEFAULT               DESCRIPTION
% d          array (up to 3d)           raw data (time down columns);
%                                       each column is treated as a
%                                       separate (time) series
% si         a) scalar                  sampling interval in µs
%            b) string 'idx'            sampling interval is set to 1; time
%                                       stamps will be given as indices
% thresh     array                      thresholds (same number of columns
%                                       and 'slices' as in d); if negative,
%                                       assume that negative-going crossing
%                                       is looked for
% detMode    char, 'cross'              - 'cross': time stamp is set to
%                                       crossing of threshold
%                                       - 'peak': each peak beyond
%                                       threshold will be assigned a time
%                                       stamp
%                  <<< OUTPUT VARIABLES <<<
%
% NAME       TYPE/DEFAULT               DESCRIPTION
% tsl        cell array of 1d arr       time stamp lists (ms or indices);
%                                       same dimension as thresh

% (c) Harald Hentschke, University of Tübingen
% Overhauled October 2017 (introduced peaks as ts)

% improvements: - allow more than one threshold per channel
detMode='cross';
pvpmod(varargin,{'detMode'});
[n1,n2,n3]=size(d);

if n2~=size(thresh,2) || n3~=size(thresh,3)
  error('thresholds must be arranged in a 1-by-N-by-M array with N,M equal to the corresponding dimensions of d');
end
if ischar(si)
  if strcmpi(si,'idx')
    si=1;
  else
    error('si must be a scalar or string ''idx''');
  end
else
  % convert to ms
  si=si*.001;
end

switch detMode
  case 'cross'
    % do it column by column (less memory consumption)
    for cnt2=1:n3
      for cnt1=1:n2
        % if thresh is negative, assume that negative-going crossing is looked for
        if thresh(1,cnt1,cnt2)<0.0
          tc=d(:,cnt1,cnt2)<=thresh(1,cnt1,cnt2);
        else
          tc=d(:,cnt1,cnt2)>=thresh(1,cnt1,cnt2);
        end
        tc=diff(tc);
        % thus, first tick above thresh is defined as time of crossing
        tsl{1,cnt1,cnt2}=(find(tc==1)+1)*si;
      end
    end
  case 'peak'
    % dd: positive-going peaks = -1, negative-going peaks = 1
    dd=diff(diff(d)>0);
    for cnt2=1:n3
      for cnt1=1:n2
        % if thresh is negative, assume that negative-going crossing is looked for
        if thresh(1,cnt1,cnt2)<0.0
          % indexes to negative-going peaks
          ix=find(dd(:,cnt1,cnt2)==1)+1;
          % pick those below threshold
          tsl{1,cnt1,cnt2}=ix(d(ix)<thresh)*si;
        else
          % indexes to positive-going peaks
          ix=find(dd(:,cnt1,cnt2)==-1)+1;
          % pick those above threshold
          tsl{1,cnt1,cnt2}=ix(d(ix)>thresh)*si;
        end
      end
    end
  otherwise
    error('bad detMode')
end

