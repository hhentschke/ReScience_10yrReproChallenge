function rmouse_specgram
% ** function rmouse_specgram
% computes & plots spectrogram of principal channel

% improvements:
% - index segments with artifacts and delete them?

global WP DS AP 


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

% scaling factor for computing psd
scalfac=AP.ppSeg/sum(window.^2);

% ------- PART II: the works
% load data
if exist([DS.dpath '\' DS.abfFn '.mat'],'file')
  d=matDload([DS.dpath '\' DS.abfFn '.mat'],'start',.001*WP.bor,'stop',.001*WP.eor,...
    'channels',AP.rawChPrincNm,'verbose',0);
elseif exist([DS.dpath '\' DS.abfFn '.raw'],'file')
  rawFInfo=rawinfo('filename',[DS.dpath '\' DS.abfFn '.raw'],'print','no');
  d=rawload([DS.dpath '\' DS.abfFn '.raw'],AP.rawChPrincNm,...
    [WP.bor  WP.eor],rawFInfo);
  % put into array and convert to mV
  d=cat(2,d{:})/1000;
  tmp=rawfi2genfi(rawFInfo,'channels',AP.rawChAnNm);
  si=tmp.si;
else
  d=abfload([DS.dpath '\' DS.abfFn '.abf'],'start',.001*WP.bor,'stop',.001*WP.eor,...
    'channels',AP.rawChPrincNm,'verbose',0);
end
if DS.rawSignalInverted
  d=-1*d;
end

% ** spectrogram.m is the newer version of specgram.m available in Matlab V
% 7.1 (?) and higher
[spg,F,T]=spectrogram(d,window,AP.dftOlapPts,AP.ppSeg,1e6/WP.osi);
% [spg,F,T]=specgram(d,AP.ppSeg,1e6/WP.osi,window,AP.dftOlapPts);

% shift time such that time points represent midpoints of intervals
T=T+discrete2cont(AP.ppSeg/2,WP.osi/1e6);
% set appropriate units for time
if diff(T([1 end]))>300
  tDiv=60;
  T=T/tDiv;
  tbase='(min)';
else
  tDiv=1;
  tbase='(s)';
end

fbinw=diff(F(1:2));
% indices into F (freq bands)
freqIx=F>=specgramFreqWin(1) & F<=specgramFreqWin(2);
% cut down and convert to power spectral density
F=F(freqIx);
spg=(spg(freqIx,:).*conj(spg(freqIx,:)))*scalfac*2/AP.ppSeg*WP.osi*1e-6;

% for sure there will be outliers (single spectra with huge power values because
% of artifacts) which mess up scaling: let's find the nasty spikes and set upper 
% color limit accordingly
co=cumh(spg(:),.0001,'p',[.01 .99]);
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

% first plot: overview
subplot(2,1,1)
imagesc(T,F,spg,clim);
set(gca,'ydir','normal');
ylabel('Frequency (Hz)');

% second plot: excerpt chosen for analysis
subplot(2,1,2)
excIx=T*tDiv>=WP.boe*.001 & T*tDiv<=WP.eoe*.001;
T=T(excIx);
spg=spg(:,excIx);
% re-compute color limit for excerpt
co=cumh(spg(:),.0001,'p',[.01 .99]);
clim=co';
imagesc(T,F,spg,clim);
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
    print(pa,[WP.figName '_spectrogram' ext]);   
  end
end

% close(fh_specg);



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
