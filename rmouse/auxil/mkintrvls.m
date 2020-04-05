function [intrvls,intrvlsPts,ilenChosen,nilenPtsChosen,remPts]=mkintrvls(intv,varargin)
% ** function [intrvls,intrvlsPts,ilenChosen,nilenPtsChosen,remPts]=mkintrvls(intv,varargin)
% divides interval intv (e.g. a period of time) into regularly spaced
% smaller sub-intervals possibly overlapping each other and puts out start
% and stop points of the intervals. For some variables, a range of
% acceptable values can be specified, and mkintrvls will select the lowest
% value resulting in the least omission of points at the right border.
%
% --------------------->>> INPUT VARIABLES >>>-----------------------------
% NAME         TYPE/DEFAULT        DESCRIPTION
% intv         2-element-array     start and end times of interval to be 
%                                  subdivided (arbitrary continuous unit, 
%                                  e.g. ms)
% resol        scalar, 1           desired resolution, e.g. sampling 
%                                  interval for time data
% ilen         scalar/2-array      length of intervals (specify fixed value
%                                  or range)
% ni           scalar/2-array      number of intervals (specify fixed value
%                                  or range)
% olap         scalar, 0           extent of overlap between segments 
% border       string, 'skip'      deal with points possibly remaining from 
%                                  division of interval: 
%                                  'skip' - ignore, all intervals will have 
%                                    identical length, but interval may not 
%                                    be completely covered
%                                  'include' - include, last interval may 
%                                    be shorter than other intervals
%                                  'embrace' - include, last interval may 
%                                    be as long as others and shifted back-
%                                    wards, resulting in more overlap be-
%                                    tween last and last but one interval
%
% -------------------- <<< OUTPUT VARIABLES <<< ---------------------------
% NAME               TYPE/DEFAULT   DESCRIPTION
% intrvls            2d-array       sampling-aligned start (col 1) and stop 
%                                   (col2) times of intervals
% intrvlsPts         2d-array       start (col 1) and stop (col2) times of 
%                                   intervals in points
% ilenChosen         scalar         chosen value of ilen, which may differ 
%                                   from input value, see above
% nilenPtsChosen     scalar         chosen value of ilenPts, which may
%                                   differ from input value, see above
% remPts             scalar         remainder (pts)


% Function renovated Nov. 2017

resol=1;
ilen=nan;
ni=nan;
olap=0;
border='skip';

pvpmod(varargin,{'resol','ilen','ni','olap','border'});

idx=1;

if numel(intv)~=2
  error('intv must be an array with two elements');
end
if any(isnan(ilen)) && any(isnan(ni))
  error('specify length of intervals or number of intervals')
end
if ~all(isnan(ilen)) && ~all(isnan(ni))
  error('specify EITHER length of intervals OR number of intervals, not both'); 
end
if ~any(isnan(ilen))
  % remember trivial but easy-to-be-unaware-of differences between continuous and discrete variables
  % - length of a time interval (continuous vars), e.g. len=diff(intrvl)
  % - length of a time interval (discrete vars, i.e. points),  len=diff(intrvlsPts)+1
  % conversion of continuous variable into discrete variable (points) is done in a 'conservative' way:
  % the discrete variable re-converted into continuous format will be <= original continuous variable
  idiff_pts=unit2pts(ilen,resol,'round')-1;
  olapdiff_pts=unit2pts(olap,resol,'round')-1;
  Int_pts(1)=unit2pts(intv(1),resol,'floor')+1;
  Int_pts(2)=unit2pts(intv(2),resol,'ceil');
  % total # of pts in intv
  IntLen_pts=diff(Int_pts)+1;
  % simple checks
  if isempty(intersect(numel(ilen),[1 2]))
    error('ilen must be a scalar or a 2-element array');
  end
  if any(idiff_pts>IntLen_pts)
    error('(range of) desired sub-interval length is longer than intv');
  end
  if any(ilen<resol)
    error('(range of) desired sub-interval length is too short given resolution'); 
  end  
  if olapdiff_pts>=idiff_pts
    error('overlap between intervals is larger than or equal to sub-interval length'); 
  end
  
  % find best interval length if a range is given
  if sum(size(idiff_pts))==3
    % idiff_pts is now an n-element 2d array
    idiff_pts=idiff_pts(1):idiff_pts(2);
  end
  [~,idx]=min(rem(IntLen_pts-olapdiff_pts,idiff_pts-olapdiff_pts));
  % now idiff_pts is a scalar whose value represents the optimal subinterval length
  idiff_pts=idiff_pts(idx);
  [remPts,~]=min(rem(IntLen_pts-olapdiff_pts-1,idiff_pts-olapdiff_pts));
  % calculate start points
  switch border
    case 'skip'
      % 'remainder is NOT included here
      intrvlsPts=[Int_pts(1):idiff_pts-olapdiff_pts:Int_pts(2)-idiff_pts]';
    case 'include'
      % remainder is included in subinterval of smaller length here
      intrvlsPts=[Int_pts(1):idiff_pts-olapdiff_pts:Int_pts(2)-olapdiff_pts-1]';
    case 'embrace'
      % remainder is included in subinterval of full length shifted such that its end falls upon end of whole data segment
      intrvlsPts=unique([Int_pts(1):idiff_pts-olapdiff_pts:Int_pts(2)-idiff_pts Int_pts(2)-idiff_pts]');
    otherwise
      error('bad input string for ''borders''');
  end
  % stop points
  intrvlsPts(:,2)=intrvlsPts(:,1)+idiff_pts;
  % for the 'include' case above check/adjust stop point of last interval
  if intrvlsPts(end,2)>Int_pts(2)
    if strcmpi(border,'include')
      intrvlsPts(end,2)=Int_pts(2);
    else
      error('internal error: last point of last interval beyond bounds of data segment');
    end
  end
  % subtract 1 to have 1st data point at continuous time 0
  intrvls=pts2unit(intrvlsPts-1,resol);   % ms
  ni=size(intrvlsPts,1);
  %   disp(['interval ' num2str(intv) ' was divided into ' int2str(ni) ' intervals of regular length '...
  %     num2str(pts2unit(diff(intrvlsPts(1,:))+1,resol)) ', remainder ' num2str(pts2unit(remPts,resol))]);
  %   disp(['in points: interval ' int2str(Int_pts) ', regular interval length '...
  %     int2str(diff(intrvlsPts(1,:))+1) ', remainder ' int2str(remPts)]);
  disp(['number of subintervals: ' int2str(ni) ', remainder: ' int2str(remPts)]);
elseif ~isnan(ni)
  error('sorry, not yet finished with this part..');
end

nilenPtsChosen=idiff_pts(1,:)+1;
ilenChosen=pts2unit(nilenPtsChosen,resol);

% ----------------- local func ------------------------------------
% assume all continuous measures have same units
function xu=pts2unit(x,resol)
xu=x*resol;

function xpts=unit2pts(x,resol,rdmeth)
tmp=x/resol;
eval(['xpts=' rdmeth '(tmp);']);
