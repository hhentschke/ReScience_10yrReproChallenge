function r=makerow(m)

% r=makerow(m)
%
% generates a row vector from matrix m by concatenating rows of m
% returns the string 'error' if ndims(m)>2

% C.Schwarz 1/2001


d=ndims(m);
if d>2
  fprintf(1,'error makecol: dimension of input variable > 2!');
  r='error';
  return
end

r=reshape(m', [1 prod(size(m))]);
