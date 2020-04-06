function [rd,varargout]=regroup(d)
% ** function [rd,varargout]=regroup(d)
% rearranges 2d-array d into nd-array rd. d contains data in a format
% suitable for use with ANOVA functions: 
% last column:       dependent variable (=observations)
% all other columns: independent variables (=levels of factors)
% ** important: the last but one column MUST be the 'subject' factor! **
% Thus, the number of columns of d corresponds to (number of factors + 1).
% The data will be rearranged in a way that allows computing means and
% standard deviations across subjects.
% Example: consider 
% d=
%          0    1.0000    1.0000    0.0027
%     0.2000    1.0000    1.0000    0.0059
%          0    1.0000    2.0000    0.0095
%     0.2000    1.0000    2.0000    0.0081
%          0    1.0000    3.0000    0.0233
%     0.2000    1.0000    3.0000    0.0127
%          0    1.0000    5.0000    0.0023
%     0.2000    1.0000    5.0000    0.0035
%          0    2.0000    1.0000    0.0015
%     0.2000    2.0000    1.0000    0.0031
%          0    2.0000    2.0000    0.0015
%     0.2000    2.0000    2.0000    0.0030
%          0    2.0000    3.0000    0.0031
%     0.2000    2.0000    3.0000    0.0037
%          0    2.0000    5.0000    0.0014
%     0.2000    2.0000    5.0000    0.0020
% 
% Let's say the data represent the following:
% first column: concentration of sugar in feeding solution
% second column: color of feeding solution (1=colorless, 2=red)
% third column: ant ID
% last column: weight (g) the respective ant can carry after being fed
% 
% The output, rd, will be a 2 by 2 by 4 array:
% rd(:,:,1) =
%     0.0027    0.0015
%     0.0059    0.0031
% rd(:,:,2) =
%     0.0095    0.0015
%     0.0081    0.0030
% rd(:,:,3) =
%     0.0233    0.0031
%     0.0127    0.0037
% rd(:,:,4) =
%     0.0023    0.0014
%     0.0035    0.0020
% 
% Thus, the number of rows, columns, etc. corresponds to the number of
% levels of the first, second, etc. factor in d. This way, we can easily
% compute means etc. of the data across subjects and plot them.
%
% The variable output arguments contain the levels of the factors in
% column arrays. The idea here is that we may want to plot the data with
% the different factors' levels on the abscissa.
% 
% Additional notes: 
% - missing data will be set to nan
% - the routine does not handle multiple observations for each unique set
%   of independent variables (repetitions of measurements)

% --- i. determine number of factors and levels & input checks
[n1 n2]=size(d);
nFactor=n2-1;
if nFactor<2
  error(['there is only one factor; running ' mfilename ' on the data is pointless']);
end

for g=1:nFactor
  % watch out - for unique nans are not unique
  fac(g).level=unique(d(~isnan(d(:,g)),g));
  fac(g).nLevel=length(fac(g).level);
  if fac(g).nLevel<2
    error(['factor # ' int2str(g) ' has only one level']);
  end
  % var output: levels
  varargout{g}=fac(g).level;
  % ** replace levels by ranks **
  for h=1:fac(g).nLevel
    d(d(:,g)==fac(g).level(h),g)=h;
  end
end

% --- ii. preallocate output array
rd=repmat(nan,[fac.nLevel]);

% --- iii. fill output arrays
% there are as many theoretical combinations as the product of all factors'
% levels
nCombin=prod([fac.nLevel]);
% string to be put into eval as output arguments for ind2sub
evalStr=[];
for g=1:nFactor
  evalStr=[evalStr 'i' int2str(g) ' '];
end
% generates variables i1, i2, ...
eval(['[' evalStr ']=ind2sub([fac.nLevel],(1:nCombin)'');']);
% ..the concatenation of which results in an array listing all possible
% combinations of the different factors' levels
eval(['combo=[' evalStr '];']);
% ..and which also represent proper indices into rd:
rdIndexStr=[int2str(combo(:,1))];
for g=2:nFactor
  rdIndexStr=[rdIndexStr  repmat(',',nCombin,1)  int2str(combo(:,g))];
end
    
% finally, fill output array
dIx=repmat(nan,nCombin,1);
for g=1:nCombin
  ix=find(~any(repmat(combo(g,:),[n1 1])-d(:,1:nFactor),2));
  if ~isempty(ix)
    if length(ix)>1
      error([mfilename ' does not handle multiple observations']);
    end
    eval(['rd(' rdIndexStr(g,:) ')=d(ix,end);']);
    dIx(g)=ix;
  end
end 