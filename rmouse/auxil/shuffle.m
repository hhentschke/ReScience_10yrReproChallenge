function [d,i]=shuffle(dim,d,varargin)
% ** function d=shuffle(dim,d,varargin)
%  shuffles elements in array d (=rearranges their positions in a pseudorandom way) 
%
%                         >>> INPUT VARIABLES >>>
%
% NAME              TYPE/DEFAULT        DESCRIPTION
% d                 array (up to 3D)    data to be shuffled
% dim               scalar              dimension along which to shuffle
%                                       (1=rows, 2=columns,...
% seed              scalar integer      seed of rand
%                     
%                         <<< OUTPUT VARIABLES <<<
%
% NAME             TYPE/DEFAULT         DESCRIPTION
% d                nd-array             shuffled data 

seed=rand('twister');
pvpmod(varargin);
rand('twister',seed);

[n1 n2 n3]=size(d);
sz=[n1 n2 n3];
dimix=setdiff(1:3,dim);

% set up strings to be eval-ed:
% (i) original index to d
dstr{dimix(1)}='ix1';
dstr{dimix(2)}='ix2';
dstr{dim}=':';
dstr=[dstr{1} ',' dstr{2} ',' dstr{3}];
% (ii) shstr contains reference to ix, the array of randomly shuffled indices
shstr{dimix(1)}='ix1';
shstr{dimix(2)}='ix2';
shstr{dim}='ix';
shstr=[shstr{1} ',' shstr{2} ',' shstr{3}];

for ix2=1:sz(dimix(2))
  for ix1=1:sz(dimix(1))
    [nada,ix]=sort(rand(size(d,dim),1));
    eval(['d(' dstr ')=d(' shstr ');']);
  end
end



% eval([' shuffd=d('  repmat([':,'],1,dim-1) 'i' repmat([',:'],1,nd-dim) ');']);