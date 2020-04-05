function d=strmread(varargin)
% ** function d=strmread(varargin)
% counterpart to strmwrite
global DS
d=i16load(varargin{:})/2^15*DS.nsRng;





