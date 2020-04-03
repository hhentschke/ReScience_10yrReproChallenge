% This is simple script pulling summary behavioral data out of a single rmouse
% data file. Each time before invoking this script, you have to run the
% anpar-generating m-file corresponding to the data file of interest. 
% ** MAKE SURE THE TIME INTERVAL (AP.rawExcerpt) IS SET CORRECTLY **

global ANPAR DSET AP DS
DSET=DS;
ANPAR=AP;
combine_bv('resol',.2);