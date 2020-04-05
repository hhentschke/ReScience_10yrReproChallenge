function [c,nc]=combin(n,varargin)
% ** function [c,nc]=combin(n,varargin) computes all non-redundant pairwise
% combinations of numbers 1:n and puts them out in a 2d-array. 
% Default order of pairs is 1-1, 1-2, ..., 1-n, 2-2,...,n-n. 
% Optional input parameters must be specified as parameter/value pairs,
% e.g. as in
%          [c,nc]=combin(n,'autoC',0);
%
%                      >>> INPUT VARIABLES >>>
% NAME          TYPE/DEFAULT          DESCRIPTION
% n             scalar                number of channels/elements
% autoC         logical, false        if true auto-pairs like [1 1] 
%                                      will be included
% diagMode      logical, false        if true order of pairs will be 
%                                      1-1, 2-2,...,n-n, 1-2, 2-3,...
%                                      (as if going along diagonals of an n
%                                      by n matrix and pairing row and
%                                      column indexes)
%                      <<< OUTPUT VARIABLES <<<
% NAME          TYPE/DEFAULT          DESCRIPTION
% c             array of size [nc 2]  the combinations
% nc            scalar                the number of combinations
%
% 
% ** THIS IS A DINOSAUR - CONSIDER USING NCHOOSEK INSTEAD **

% default values
autoC=false;
diagMode=false;
pvpmod(varargin);

if autoC
  nc=sum(1:n);
  offset=0;
else
  nc=sum(1:n)-n;
  offset=1;
end

if diagMode
  c=[];
  m=reshape(1:n^2,n,n);
  for dIx=1+offset:n
    [rIx,cIx]=ind2sub(n,diag(m,dIx-1));
    c=cat(1,c,[rIx cIx]);
  end
else
  c=zeros(nc,2);
  idx=0;
  for i=1:n
    for g=i+offset:n
      idx=idx+1;
      c(idx,:)=[i g];
    end
  end
end

  