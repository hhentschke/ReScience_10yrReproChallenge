function rcatfish(rv,varargin)
% ** function rcatfish(rv,varargin)
% fishes out a subset of concatenated channel pair data from struct R as produced by
% rcat.m, brings them in a format suitable for principal components analysis and puts
% them in global variable D
% ** expects global var R in base workspace **
%                         >>> INPUT VARIABLES >>>
% NAME           TYPE/DEFAULT          DESCRIPTION
% rv             cell array of chars   any (combination) of fields of R containing
%                                      segment-wise computed parameters, like 'thCCPeak'
% chanComb       char arr, 'neigh'     'neigh' - nearest neighbor (line plots)
%                                      'princ' - all vs. principal channel (line plots)
%                                      'all'   - all, color-coded matrix plot
%                                      'coef'  - only those not crosscorrelating with
%                                                each other significantly (NYI)
% memEcon        scalar, 0             if nonzero, fields of R will be deleted as soon as
%                                      they are not needed anymore in this function (because 
%                                      both old fields of R and the new field may take up
%                                      lots of memory

global DS AP R D

% to do: 
% - implement choice based on corrcoef
% - implement memEcon
% - ignore chanComb if thgaeCC ismember of rv
% - option to restrict range of channels but keeping the architecture of channel combos
% - spit out meaningful names (e.g. thCCPeak_100p) for each column of the new field so
%   they can be easily identified


% --- defaults
pvpmod(varargin);

% --- preliminaries: local copies of AP and DS, checks & some vars
AP=R.AP;
DS=R.DS;
rmouse_APcheck;
rawCh=rmouse_chan;
nrv=length(rv);
nt=size(R.tb,1);

% @@ thgaeCC?
switch chanComb
  case 'neigh'
    template1=repmat(nan,nt,nLFPCh-1);    
  case 'princ'
    template1=repmat(nan,nt,nLFPCh);    
  case 'all'
    template1=repmat(nan,nt,(nLFPCh-1)*nLFPCh/2);
  case 'coef'    
    error('''coef'' not yet implemented');    
  otherwise
    error('bad choice of chanComb');
end


% --- check whether all fields requested exist
nixex=[]; nixi=[];
fina=fieldnames(R);
for rvi=1:nrv
  if isempty(strmatch(rv{rvi},fina))
    nixi=[nixi rvi];
    nixex=[nixex ' ' rv{rvi}];
  end
end
if length(nixi) ~= 0
  error(['sorry, pal, field(s) ' nixex ' do not exist in R']);
end
clear nix*

% --- reshape matrices as linear (1D) arrays, time running down columns & combine 
disp('** fishing for & combining columns..');
D=[];
for rvi=1:nrv
  y=template1;
  % select extraction method depending on results var chosen (in essence,
  % differentiate between 1 by n by t and n by n by t
  % also, make sure autocorrs are not included
  switch rv{rvi}
    case {'thgaeCCPeak','thgaeCCPeakT'}
    otherwise
      switch chanComb
        case 'neigh'
          ct=0;
          for ct=1:nLFPCh-1
            coli=AP.LFPccInd(ct)+1;
            eval(['y(:,ct)=permute(R.' rv{rvi} '(AP.LFPccInd(ct),coli,:),[3 2 1]);']);
          end
        case 'princ'
          for ct=1:nLFPCh
            eval(['y(:,ct)=permute(R.' rv{rvi} '(min(AP.LFPccInd(ct),AP.LFPpcInd2),max(AP.LFPccInd(ct),AP.LFPpcInd2),:),[3 2 1]);']);
          end
        case 'all'
          ct=0;
%           for rowi=1:nLFPCh-1
%             for coli=rowi+1:nLFPCh
          for rowi=AP.LFPccInd(1:end-1)
            for coli=AP.LFPccInd(2:end)
              ct=ct+1;
              eval(['y(:,ct)=permute(R.' rv{rvi} '(rowi,coli,:),[3 2 1]);']);
            end
          end
        case 'coef'    
      end
  end % switch rv{rvi}
  D=cat(2,D,y);
end % for rvi=1:nrv

clear y template1
