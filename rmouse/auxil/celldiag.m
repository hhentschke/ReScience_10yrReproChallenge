function x = celldiag(d, varargin)
% x = celldiag(varargin)
% Wrapper for function diag for cell arrays
[r, c] = size(d);
linix_array = reshape(1:r*c, [r, c]);
if nargin == 1
    ix = diag(linix_array);
else
    ix = diag(linix_array, varargin{1});
end
x = d(ix);