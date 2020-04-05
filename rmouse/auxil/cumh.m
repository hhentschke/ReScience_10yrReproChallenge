function [co,ncadh,bins]=cumh(d,res,varargin)
% **function [co,ncadh,bins]=cumh(d,res,varargin)
%  produces a normalized cumulative amplitude histogram of data d plus derived
%  measures (see below)
%                         >>> INPUT VARIABLES >>>
%
% NAME        TYPE/DEFAULT          DESCRIPTION
% d           array (up to 3D)      data, arranged in columns
% res         scalar                resolution of the histogram, expressed as fraction 
%                                    of the amplitude range defined by absolute  
%                                    extrema in d (e.g., .01 means 1 % resolution)
%                                    ** NOTE: if columns diverge widely in their
%                                    amplitude maxima run this routine on each 
%                                    column separately
% p           1d-array,[.01 .99]    percentiles that shall be computed (0<p<1, so
%                                    .01 is the first percentile)
%                                            
%                         <<< OUTPUT VARIABLES <<<
%
% NAME        TYPE/DEFAULT           DESCRIPTION
% co          array                  'cutoff values': the percentiles 
%                                      corresponding to values given in p
%                                      (in columns corresponding to those in d)
% ncadh       array (ndim as d)       normalized (to 1) cumulative amplitude
%                                      distribution histogram of d
% bins        array                   centers of bins

warning('THIS FUNCTION IS A DINOSAUR - consider using more adequate functions (e.g. prctile)')

% ----- default values & varargin -----
verbose=0;
p=[.01 .99];
pvpmod(varargin);
if verbose, disp(['** ' mfilename]); end;
if res<=0 | res>=1, error('check resol'); end
if any(p<0 | p>1), error('check p'); end

[n1,n2,n3]=size(d);
% get extreme values
ma=max(max(max(d,[],1),[],2),[],3);
mi=min(min(min(d,[],1),[],2),[],3);
if ma==mi
  warning('elements of d are identical - output all set to []');
  co=[];
  ncadh=[];
  bins=[];
else
  % resolution: make such that last 'edge' value will be 1/1000 larger
  % than largest value in d (so that it will not reside lonely in last bin)
  resol=(ma-mi)*res*1.001;
  % edges array
  edges=mi:resol:ma+resol;
  % centers of bins
  bins=edges(1:end-1)+resol*.05;
  
  n=histc(d,edges,1);
  % remove last row of n (it is zero)
  n(end,:,:)=[];
  ncadh=cumsum(n);
  % normalized cumulative amplitude distribution
  ncadh=ncadh./repmat(ncadh(end,:,:),[length(bins) 1 1]);
  
  
  for i=1:length(p)
    [nix,ix]=min(abs(ncadh-p(i)),[],1);
    co(i,1:n2,1:n3)=bins(ix);
  end
end