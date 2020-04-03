function R=combine_rr02
% extracts and combines 'raw results' (rr) like spectra or crosscorrelations 
% from a collection of data sets (experiments).
% needs ANPAR and DSET (=concatenation of individual and matching AP and DS)
% one experiment per column. 
%
% combine_rr works similar to the combineX (X any number from 1 to 4) functions
% in that it loops over data sets and collects data. The major difference is the nature of
% the data collected (single parameters in the other combine funcs, whole series here) and
% as a consequence, the structure of the final results variable R

% to do: get rid of 60 Hz!

global ANPAR DSET

% I. Set parameters ---------------------------------------------------------
% choose single results variable to collect and average/plot - 
% must be a field name of r (will be put in eval)
% currently, works only with power spectra (rawPMn)
rv={'rawPMn'};
% choose auto- or cross-channel results (cross not computed for theta, gamma env corr)
q='cross';
q='auto';
% choose behaviors to be treated (legal value of AP.segmentType)
behav={'immobile','exploring'};
nBehav=length(behav);

pcol={'b','r'};

% spectra etc. will be up/downsampled to this resolution (Hz)
% -> needed because sampling intervals were not identical in all cases
fResol=.25;
% desired frequency range of spectra 
rfreq=[40:fResol:90];
% rfreq=[5:fResol:30];
% frequency range which shall be considered for normalization (must be contained in above)
normrfreq=rfreq; 
% kind of normalization of (difference) spectra:
% a) 'none'
% b) 'ctrlPow': for each channel, divide spectra by total spectral power in control case 
%    (= divide by integral (=area) of control spectrum)
%    -> this normalization emphasizes RELATIVE differences between control and drug
%    condition, even if the absolute differences are small
% c) 'diffPow': for each channel, normalize DIFFERENCE of spectra such that their area is equally 1
%    -> this normalization upscales channels with low absolute power
normType='ctrlPow';
normType='diffPow';
normType='none';

% shall log of power density be plotted?
logp=0;
% print?
printas=[];'-djpeg90';

% II. Various preparations ---------------------------------------------------------
rfreqix=1:length(rfreq);
normrfreqix=find(rfreq>=normrfreq(1) & rfreq<=normrfreq(end));

rmouse_ini;

% --- plot stuff
% this generates variable segTypeGlobP containing global settings for plots of behaviors
stgp;
% in plots, set vertical lines at these frequency values
vlines=[6 8 20 40];
close all;
labelscale('fontSz',10,'scaleFac',1.0,'lineW',1.0,'markSz',6); 

% -------- PART III: collect information about recording sites & set up results var
[n1 n2 n3]=size(ANPAR);
[dn1 dn2 dn3]=size(DSET);
if ~isequal([n1 n2 n3],[dn1 dn2 dn3]), error('DSET and ANPAR must have equal dimension'); end


% variable holding collected results: R
% - a struct array (as many elements as genotypes)
%   - freq down columns, 
%   - sessions along row
%   - different behaviors in slices
% freq is unknown at this point
for i=1:n3
  R(i).r=repmat(nan,[1 n2 n3]);
end

% -------- PART IV: collection of data
% loop over data sets: 
ri=1;
for i3=1:n3
  for ci=1:n2
    AP=ANPAR(ri,ci,i3);
    DS=DSET(ri,ci,i3);
    % need rmouse_chan only for WP.elx
    rawCh=rmouse_chan;
    % load results var..
    if isempty(strfind(AP.resPath,':')), AP.resPath=[WP.rootPath AP.resPath]; end    
    load([AP.resPath '\' AP.resFn],'r');
    % ..find behaviors..
    for bi=1:length(behav)
      % bix is index into r
      bix(bi)=strmatch(behav{bi},{r(:).segmentType});
    end
    % find freqs of interest (end values must be a little wider because interpolation
    % will yield nans for x values outside interval)
    frix=find(r(1).F>=rfreq(1)-fResol & r(1).F<=rfreq(end)+fResol);
    % preparations for cutting out line hum (from complete spectrum)
    humF=60:60:r(1).F(end);
    for h=1:length(humF)
      % indices to immediately adjacent bins
      [nada,tmpix]=min(abs(r(1).F-humF(h)));
      % the ones to be replaced
      ix1a=tmpix-1:tmpix+1;
      ix1a(ix1a<1 | ix1a>length(r(1).F))=[];
      ix1{h}=ix1a;
      % the ones to compute mean from & replace with 
      ix2a=[tmpix-5:tmpix-2 tmpix+2:tmpix+5];
      ix2a(ix2a<1 | ix2a>length(r(1).F))=[];            
      ix2{h}=ix2a;
    end
    % and extract data
    rvi=1;
    for bi=1:length(bix)
      % tell what we're dealing with
      disp(['data: ' AP.resFn ', behavior: ' behav{bi} ', data:' rv{rvi}]);
      % tmpr is the original data
      eval(['tmpr=r(bix(bi)).' rv{rvi} ';']);
      allEmpty=isequal(size(tmpr),[1 1]);
      % extract data, interpolate and put in proper position
      if allEmpty
        % this accounts for missing behaviors: fill with nans
        disp('**** no data, filling with nans');
        R(i3).r(rfreqix,ci,bi)=repmat(nan,rfreqix(end),1);
      else
        % whole spectrum of current channel pair
        % ** note: for cross measures, compute absolute value here
        tmpP=abs(tmpr{AP.LFPpcInd2,AP.LFPpcInd2});
        for h=1:length(ix1)
          % this substitutes 60 Hz peaks and harmonics by mean of neighboring bins
          tmpP(ix1{h},:)=repmat(mean(tmpP(ix2{h},:),1),[3 1]);
        end
        R(i3).r(rfreqix,ci,bi)=interp1(r(1).F(frix),tmpP(frix),rfreq);
      end
    end % for:behav
  end % for:cols of ANPAR=experiments
end % for:slices of ANPAR=genotypes




fh=mkfig('huhu');
for i3=1:n3
  cuh=cumsum(R(i3).r);
  cuh=cuh./repmat(cuh(end,:,:),[length(rfreq) 1 1]);
  for bi=1:length(bix)
    subplot(2,2,bi);
    hold on
    plot(rfreq,R(i3).r(:,:,bi),pcol{i3});
    subplot(2,2,2+bi);
    hold on
    plot(rfreq,cuh(:,:,bi),pcol{i3});
  end
end
    

% ----- local funcs -----------

function figha=mkfig(ftag)
figha=findobj('tag',ftag);
if isempty(figha), figha=figure;
else  figure(figha);
end
tmpScrSz=get(0,'Screensize')*.8+.1*rand;
tmpScrSz([1 2])=tmpScrSz([1 2])+tmpScrSz([3 4])*.10;
set(figha,'position',tmpScrSz,'tag',ftag,'name',ftag,...
  'color',[.9 .9 1],'numbertitle','off');
clf;
