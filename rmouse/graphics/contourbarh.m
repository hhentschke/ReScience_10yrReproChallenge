function lh=contourBarh(x,y,varargin)
% ** function ph=contourBarh(x,y,varargin)
% draws a HORIZONTAL bar contour of y vs x and returns handle to line.
% varargin is the string of properties passed to 'line'.
% Accepts only regularly spaced data points (as produced by e.g. hist).
% The present version may be the basis of a more elaborate function as an
% alternative to the Matlab stairs function.

tmp=unique(diff(x));
tmpmax=max(tmp); tmpmin=min(tmp);
if tmpmax/tmpmin >1.01,
   warning('x is not strictly monotonously increasing');
end;

hdist=(x(2)-x(1))*.5;
xx(1:2:2*length(x)-1)=x-hdist;
xx(2:2:2*length(x))=x+hdist;
yy(1:2:2*length(x)-1)=y;
yy(2:2:2*length(x))=y;

lh=line(yy,xx);
set(lh,varargin{:});
