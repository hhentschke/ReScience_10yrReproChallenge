% computes & plots spectrogram of principal channel

% improvements:
% - index segments with artifacts and delete them?

% ------ PART I: preparations
% check whether AP.ppSeg is really a power of 2 - this is not necessary here,
% but in rmouse_cspecp
if log2(AP.ppSeg)~=round(log2(AP.ppSeg)), error('AP.ppSeg must be a power of 2'); end
% window
if strcmpi(AP.dftWin,'rect'), window=ones(AP.ppSeg,1); 
else eval(['window=' AP.dftWin '(AP.ppSeg);']);
end

% frequency range to plot (should be a parameter for ANPAR)
specgramFreqWin=[1 30];

% the power spectra will be smoothed by convoluting them with a triangular window.
% the parameter below specifies the half width of that window in Hz
specGramKrnlHWid=                 .24;

% ------- PART II: the works
% load data
if exist([dpath '\' abfFn '.mat'],'file')
  d=matDload([dpath '\' abfFn '.mat'],'start',.001*boe,'stop',.001*eoe,...
    'channels',AP.rawChPrincNm,'verbose',verbose);
else
  d=abfload([dpath '\' abfFn '.abf'],'start',.001*boe,'stop',.001*eoe,...
    'channels',AP.rawChPrincNm,'verbose',verbose);
end
if DS.rawSignalInverted
  d=-1*d;
end
% Each column of B contains an estimate of the short-term, time-localized
% frequency content of the signal A.  Time increases linearly across the
% columns of B, from left to right.  Frequency increases linearly down
% the rows, starting at 0.
[spg,F,T]=specgram(d,AP.ppSeg,1e6/osi,window,AP.dftOlapPts);
% shift time such that time points represent midpoints of intervals
T=T+discrete2cont(AP.ppSeg/2,osi/1e6);
% add offset and set appropriate units for time
T=T+boe*.001;
if diff(T([1 end]))>300
  T=T/60;
  tbase='(min)';
else
  tbase='(s)';
end

fbinw=diff(F(1:2));
% indices into F (freq bands)
freqIx=F>=specgramFreqWin(1) & F<=specgramFreqWin(2);
% cut down and convert to power spectral density
F=F(freqIx);
spg=abs(spg(freqIx,:))/fbinw;

% for sure there will be outliers (single spectra with huge power values because
% of artifacts) which mess up scaling: let's find the nasty spikes and set upper 
% color limit accordingly
co=cumh(spg(:),.0001,'p',[.005 .995]);
clim=co';

% smooth (convolution with triangular kernel
% half width, in bins, of triangular kernel for smoothing of spectra 
% (more exactly, kernel will have length 2*krnlhPts+1)
krnlhPts=ceil(specGramKrnlHWid/fbinw);
krnl=triang(2*krnlhPts+1);
% don't forget to normalize
krnl=krnl/sum(krnl);
for ii=1:length(T)
  tmpSpg=conv(spg(:,ii),krnl);
  % get rid of krnlhPts border points on either end
  spg(:,ii)=tmpSpg(krnlhPts+1:end-krnlhPts);
end

% --- figure 
tmpftag='Spectrogram';
fh_specg=findobj('tag',tmpftag);
if isempty(fh_specg), fh_specg=figure;
else  figure(fh_specg);
end
tmpScrSz=get(0,'Screensize');
tmpScrSz([3 4])=round([tmpScrSz(3)*.35 tmpScrSz(4)*.8]);  
set(fh_specg,'position',tmpScrSz,'tag',tmpftag,'name','spectrogram' ,...
  'color',[0.8 0.8 0.8],'numbertitle','off');

% colormap
colormap copper;

% plot: linear scale
subplot(2,1,1)
imagesc(T,F,spg,clim);
set(gca,'ydir','normal');
ylabel('Frequency (Hz)');

% plot: log scale
subplot(2,1,2)
imagesc(T,F,log(spg));
set(gca,'ydir','normal');
ylabel('Frequency (Hz)');
xlabel(['Time ' tbase]);

% print this summary figure?
if ~isempty(AP.printas{1}), 
  for i=1:length(AP.printas)
    pa=AP.printas{i};
    if strfind(pa,'ps'), ext='.ps';
    elseif strfind(pa,'jpeg'), ext='.jpg';
    else ext='';
    end
    print(pa,[figName '_spectrogram' ext]);   
  end
end

% close(fh_specg);

% cleanup
clear tmp* window specgramFreqWin logp d spg F T tbase freqIx co clim tmpftag krnl*


% % preparations for cutting out line hum
% humF=60:60:F(end);
% for h=1:length(humF)
%   % indices to immediately adjacent bins
%   [nix,ix]=min(abs(F-humF(h)));
%   % the ones to be replaced
%   ix1a=ix-1:ix+1;
%   ix1a(ix1a<1 | ix1a>length(F))=[];
%   ix1{h}=ix1a;
%   % the ones to compute mean from & replace with 
%   ix2a=[ix-5:ix-2 ix+2:ix+5];
%   ix2a(ix2a<1 | ix2a>length(F))=[];            
%   ix2{h}=ix2a;
% end
