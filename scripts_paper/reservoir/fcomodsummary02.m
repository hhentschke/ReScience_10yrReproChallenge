function fComodSummary02
% this is a pathetic piece of code, exporting integral of pixels 
% in comodulograms corresponding to freq combinations listed in 
% combine_fcomod in a format compatible with depthp02/03

rmouse_ini;

behav={'immobile','exploring'};
compareFactor='drug';

% results variable R:
% - a struct array (as many elements as behaviors)
% - the field(s) are 1D cell arrays (one cell per channel); 
%   each data set will be embedded in the correct position
% - each cell is a 3D array:
%   - freq combo down columns, 
%   - column order control | drug | recovery (optional)
%   - different experiments in slices

% --- load 
load fComod_all;

% standard set of 7 electrodes
% ORDER MATTERS
curChInd=hp.princChInd:-1:hp.princChInd-6;

nFBand=size(fBand,1);
[combo,ncombin]=combin(nFBand,'autoC',1);

% **** make sure this is the right freq combo! ****
% theta hi-gamma
fbIx=5;
% theta lo-gamma
fbIx=3;


cnt=0;

% --- 
for bi=1:2
  kittycat=[];
  indv=[];
  for g=curChInd
    % rec depth in mm
    rs=repmat(round(g-hp.princChInd)*-.1,length(R(bi).indv{g}),1);
    % column order: rec site | control | drug
    kittycat=cat(1,kittycat, [rs   permute(R(bi).r{g}(fbIx,:,:),[3 2 1])]);
    indv=cat(1,indv, permute(R(bi).indv{g},[3 2 1]));

  end
  % create ds1, ds2, indv1, indv2
  ds1=kittycat(:,[1 2]);
  ds2=kittycat(:,[1 3]);
  indv1=indv;
  indv2=indv;
  mds1=[];
  mds2=[];
  ds1fit=[];
  ds2fit=[];
  fitx=[];
  fn=[behav{bi} '_thGaComod_auto' ];
  save([WP.rootPath '\beta3_wtko\export\' fn],'compareFactor','ds1','ds2','mds*','ds1fit','ds2fit','fitx','indv*');
  
end % for:behaviors

