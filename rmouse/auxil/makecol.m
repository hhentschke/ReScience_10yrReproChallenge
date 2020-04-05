function r=makecol(m)

% r=makecol(m)
%
% generates a col vector from matrix m by concatenating columns of m
% returns the string 'error' if ndims(m)>2

% C.Schwarz 1/2001

% d=ndims(m);
% if d>2
%   fprintf(1,'error makecol: dimension of input variable > 2!');
%   r='error';
%   return
% end
% 
% r=reshape(m, [prod(size(m)),1]);


% HH 2007
d=ndims(m);
if d>2
  error('input array must be 1D or 2D');
end
r=m(:);