function rmouse_cspecp(rawCh)
% computes power and cross power spectral density
% to do:
% - update computation of peak memory demand

global DS AP WP D r logstr


% ------ PART I: preparations
rmouse_rcleanup;

% --- check whether AP.ppSeg is really a power of 2
if log2(AP.ppSeg)~=round(log2(AP.ppSeg)), error('AP.ppSeg must be a power of 2'); end

% --- window
if strcmpi(AP.dftWin,'rect'), window=ones(AP.ppSeg,1); 
else eval(['window=' AP.dftWin '(AP.ppSeg);']);
end
% scale: sum of squared elements of window  divided by (scalfac*AP.ppSeg) 
% must be==1
scalfac=AP.ppSeg/sum(window.^2);

% --- frequency: only in 1st element of r (is same for all segmentTypes)
F=makecol(1e6/WP.osi/AP.ppSeg*[0:AP.ppSeg/2-1]);
r(1).F=F;
fbinw=diff(F(1:2));

% --- kernel for smoothing of spectra
% half width, in bins, of triangular kernel (more exactly, kernel will 
% have length 2*krnlhPts+1)
krnlhPts=ceil(AP.specKrnlHWid/fbinw);
krnl=triang(2*krnlhPts+1);
% dont forget to normalize
krnl=krnl/sum(krnl);

% --- half width of interval for computation of peak integral in bins
inthPts=round(AP.peakHW/fbinw);

% --- indices into F (freq bands) - narrow theta range will be determined
% dynamically below
deFIx=F>=AP.delta(1) & F<=AP.delta(2);
thFIx=F>=AP.theta(1) & F<=AP.theta(2);
beFIx=F>=AP.beta(1) & F<=AP.beta(2);
gaFIx=F>=AP.gamma(1) & F<=AP.gamma(2);
gaNarrowFIx=F>=AP.gammaNarrow(1) & F<=AP.gammaNarrow(2);
riFIx=F>=AP.ripple(1) & F<=AP.ripple(2);      

% --- comodulogram: 
% freq resolution shall be approximately 0.5 Hz (or larger,
% if that is the elementary freq resolution), so determine how many bins of
% the original resolution are to be summed up to to yield this
comNFBbin=max(round(0.5/fbinw),fbinw);
% freq range indices for comodulogram
tmpix=find(F>=1 & F<=100);
% starting indices
comFStartIx=max(1,tmpix(1)-floor(comNFBbin/2)):comNFBbin:min(length(F),tmpix(end)-ceil(comNFBbin/2));
% frequency bins in Hz for comodulogram (centers)
r(1).comF=F(comFStartIx)+comNFBbin/2*fbinw;

% --- index to range of freqs to plot 
prngix=find(F>=min(AP.peakFRng(1),AP.FIRcf(1)) & F<=max(AP.peakFRng(2),AP.FIRcf(2)));      

% --- peak detection: indices to spectra within which to search 
% if interval borders are too close to lower and/or upper range of freqs coerce
% to values permitting convolution AND computation of peak integral without 
% border effects 
peakix(1)=min(find(F>=AP.peakFRng(1)));
peakix(2)=max(find(F<=AP.peakFRng(2)));      
tmpflag=0;
% convolution extends length by 2*krnlhPts; however, the trailing and last
% 2*krnlhPts points (=4*krnlhPts all in all) are border points to be cut out
if peakix(1)<krnlhPts+inthPts, peakix(1)=krnlhPts+inthPts; tmpflag=1; end
if peakix(2)>length(F)-(krnlhPts-1)-inthPts, peakix(2)=length(F)-(krnlhPts-1)-inthPts; tmpflag=1; end
if tmpflag
  logstr{end+1}=['interval for peak detection in power spectra needed to be adjusted to [' num2str(makerow(F(peakix))) '] Hz'];
  warning(logstr{end});
end
% the same, but including borders
pPeakix=peakix+[-1 1]*(krnlhPts+inthPts);
% expand
pPeakix=pPeakix(1):pPeakix(2);      

% --- peak detection: DFT of the FIR filter with which to multiply averaged
%     DFT of data
% length=order of filter must be < AP.ppSeg (otherwise its DFT will be computed
% from truncated version) and should in general be fixed
bft=fft(fir1(AP.FIRfo,AP.FIRcf*1e-6*WP.osi*2),AP.ppSeg);
% keep one-sided DFT only and make a column vector
bft=(bft(1:AP.ppSeg/2))';
% 'unscaled' power is also needed:
bftP=bft.*conj(bft);

% --- preparations for cutting out line hum
humF=60:60:F(end);
for h=1:length(humF)
  % indices to immediately adjacent bins
  [foo,ix]=min(abs(F-humF(h)));
  % the ones to be replaced
  ix1a=ix-1:ix+1;
  ix1a(ix1a<1 | ix1a>length(F))=[];
  ix1{h}=ix1a;
  % the ones to compute mean from & replace with 
  ix2a=[ix-5:ix-2 ix+2:ix+5];
  ix2a(ix2a<1 | ix2a>length(F))=[];            
  ix2{h}=ix2a;
end

% --- figure
tmpftag='SpecPeakDet';
fh_spDet=findobj('tag',tmpftag);
if isempty(fh_spDet), fh_spDet=figure;
else  figure(fh_spDet);
end
tmpScrSz=get(0,'Screensize');
tmpScrSz([3 4])=round([tmpScrSz(3)*.35 tmpScrSz(4)*.8]);  
set(fh_spDet,'position',tmpScrSz,'tag',tmpftag,'name','peak detection in power spectra' ,...
  'color',[0.8 0.8 0.8],'numbertitle','off','menubar','none');
% --- plot 1: magnitude (=frequency response) of FIR filter
subplot(2,1,1)
% luckily, F computed above does the trick
plot(F(prngix),abs(bft(prngix)),'k');
niceyax;
title('frequency response FIR filter');
xlabel('freq (Hz)');
ylabel('magnitude');      

% ------------ PART II: the works
for i=1:length(r)
  if ~isempty(r(i).iPts)
    disp([r(i).segmentType '..']);
    % preallocation for comodulogram: 
    % - temporary var for reduced spectrum (order: segment | freq)
    tmpRedP=repmat(nan,r(i).ni,length(comFStartIx));
    % - correlation values:
    % frequencies in rows and columns, channels along third dim
    r(i).fComod=repmat(nan,[length(comFStartIx)*[1 1] AP.nLFPCh]);
    % - P values
    r(i).fComodP=r(i).fComod;
    % *** note that the two variables above are the only ones of the spectral analysis
    % performed here which 
    % (i) contain data only from autospectra
    % (ii) contain data 'slots' (i.e. slices) ONLY FOR ANALYZED CHANNELS, whereas in
    % all other variables non-analyzed channels are represented by nans as
    % placeholders

    % 'vector of nans'
    tmpvon=repmat(nan,1,r(i).ni);
    % template for all running pars..
    ccRunTemplate=WP.ccDerTemplate;
    % ..initialized with vector of nans
    ccRunTemplate(AP.trixie)={tmpvon};
    % assign preallocated var templates:
    % cross and power spectral density (whole freq range): averages
    r(i).rawPMn=WP.psTemplate;
    r(i).rawPStd=WP.psTemplate;
    % power in freq bands: running & averages 
    r(i).rawDePE=ccRunTemplate;
    r(i).rawDePEMn=WP.ccDerTemplate;
    r(i).rawDePEStd=WP.ccDerTemplate;
    r(i).rawThPE=ccRunTemplate;
    r(i).rawThPEMn=WP.ccDerTemplate;
    r(i).rawThPEStd=WP.ccDerTemplate;
    
    r(i).rawThNarrowPE=ccRunTemplate;
    r(i).rawThNarrowPEMn=WP.ccDerTemplate;
    r(i).rawThNarrowPEStd=WP.ccDerTemplate;

    r(i).rawBePE=ccRunTemplate;
    r(i).rawBePEMn=WP.ccDerTemplate;
    r(i).rawBePEStd=WP.ccDerTemplate;
    r(i).rawGaPE=ccRunTemplate;
    r(i).rawGaPEMn=WP.ccDerTemplate;
    r(i).rawGaPEStd=WP.ccDerTemplate;

    r(i).rawGaNarrowPE=ccRunTemplate;
    r(i).rawGaNarrowPEMn=WP.ccDerTemplate;
    r(i).rawGaNarrowPEStd=WP.ccDerTemplate;

    r(i).rawRiPE=ccRunTemplate;
    r(i).rawRiPEMn=WP.ccDerTemplate;
    r(i).rawRiPEStd=WP.ccDerTemplate;
    
%     % power in freq bands: computed from averages 
%     r(i).rawPMnDeP=WP.ccDerTemplate;
%     r(i).rawPMnThP=WP.ccDerTemplate;
%     r(i).rawPMnThNarrowP=WP.ccDerTemplate;
%     r(i).rawPMnBeP=WP.ccDerTemplate;
%     r(i).rawPMnGaP=WP.ccDerTemplate;
%     r(i).rawPMnRiP=WP.ccDerTemplate;
    
    % coherence (whole freq range): averages
    r(i).rawCohMn=WP.psTemplate;
    r(i).rawCohStd=WP.psTemplate;    
    % average coherence in frequency bands: derived from averages (see bottom of mfile)
    r(i).rawCohMnDe=WP.ccDerTemplate;
    r(i).rawCohStdDe=WP.ccDerTemplate;    
    r(i).rawCohMnTh=WP.ccDerTemplate;
    r(i).rawCohStdTh=WP.ccDerTemplate;    
    r(i).rawCohMnThNarrow=WP.ccDerTemplate;
    r(i).rawCohStdThNarrow=WP.ccDerTemplate;    
    r(i).rawCohMnBe=WP.ccDerTemplate;
    r(i).rawCohStdBe=WP.ccDerTemplate;    
    r(i).rawCohMnGa=WP.ccDerTemplate;    
    r(i).rawCohStdGa=WP.ccDerTemplate;    
    r(i).rawCohMnGaNarrow=WP.ccDerTemplate;    
    r(i).rawCohStdGaNarrow=WP.ccDerTemplate;    
    r(i).rawCohMnRi=WP.ccDerTemplate;
    r(i).rawCohStdRi=WP.ccDerTemplate;    
    % theta peak magnitude and freq: derived from average
    r(i).rawPMnPeak=WP.ccDerTemplate;
    r(i).rawPMnPeakT=WP.ccDerTemplate;    
    % theta peak magnitude and freq: running & averages
    r(i).rawPPeak=ccRunTemplate;
    r(i).rawPPeakT=ccRunTemplate;
    r(i).rawPPeakMn=WP.ccDerTemplate;
    r(i).rawPPeakStd=WP.ccDerTemplate;
    r(i).rawPPeakTMn=WP.ccDerTemplate;    
    r(i).rawPPeakTStd=WP.ccDerTemplate;    
    % spectral edge freq: running & averages
    r(i).rawPSEF=ccRunTemplate;
    r(i).rawPSEFMn=WP.ccDerTemplate;
    r(i).rawPSEFStd=WP.ccDerTemplate;
    % gamma centroid: running & averages
    r(i).rawGaCentroid=ccRunTemplate;
    r(i).rawGaCentroidMn=WP.ccDerTemplate;
    r(i).rawGaCentroidStd=WP.ccDerTemplate;
    
    % complex intermediate results var holding DFT from all LFP channels to be analyzed
    tmpDFT=repmat(nan+sqrt(-1),[AP.ppSeg/2 r(i).ni AP.nLFPCh]);
    % disp(['intermediate results var takes up ' num2str(prod(size(tmpDFT))*8/2^20,'%4.1f') ' Mb']);            
    % ------ PART IIa: DFT of each LFP channel
    for chInd=1:AP.nLFPCh
      disp([r(i).segmentType ': DFT raw ch ' rawCh(AP.LFPInd(chInd)).nm ' ..']);
      % load data
      if exist([DS.dpath '\' DS.abfFn '.mat'],'file')
        tmpd=matDload([DS.dpath '\' DS.abfFn '.mat'],'start',0,'stop',discrete2cont(r(i).iPts(end,2),WP.osi*1e-6,'intv',1),...
          'channels',{rawCh(AP.LFPInd(chInd)).nm},'verbose',WP.verbose);
      elseif exist([DS.dpath '\' DS.abfFn '.raw'],'file')
        rawFInfo=rawinfo('filename',[DS.dpath '\' DS.abfFn '.raw'],'print','no');
        tmpd=rawload([DS.dpath '\' DS.abfFn '.raw'],{rawCh(AP.LFPInd(chInd)).nm},...
          [0 discrete2cont(r(i).iPts(end,2),WP.osi*1e-3,'intv',1)],rawFInfo);
        % put into array and convert to mV
        tmpd=cat(2,tmpd{:})/1000;
        tmp=rawfi2genfi(rawFInfo,'channels',AP.rawChAnNm);
        si=tmp.si;
      else
        tmpd=abfload([DS.dpath '\' DS.abfFn '.abf'],'start',0,'stop',discrete2cont(r(i).iPts(end,2),WP.osi*1e-6,'intv',1),...
          'channels',{rawCh(AP.LFPInd(chInd)).nm},'verbose',WP.verbose);
      end
      if DS.rawSignalInverted
        tmpd=-1*tmpd;
      end
      % Due to the fact that abfload computes discrete time from continuous time a little
      % differently than is done in rmouse_xx, in some cases a single data point is
      % amiss. This is bad, and calls for a revision of abfload. In the
      % meantime, here's the workaround:
      tmppd=r(i).iPts(end,2)-length(tmpd);
      if tmppd,
        if tmppd==1
          disp('missed a point');
          tmpd(end+1)=tmpd(end);
        elseif tmppd== -1
          disp('one point too much');
        else
          error('abfload, matDload or rawload handles time information poorly');
        end
      end
      for k=r(i).ni:-1:1
        tmpIdx=r(i).iPts(k,1):r(i).iPts(k,2);
        % data segments are detrended individually
        tmpfft=fft(detrend(tmpd(tmpIdx),'constant').*window);
        % keep one-sided
        tmpDFT(:,k,chInd)=tmpfft(1:AP.ppSeg/2);
      end
    end
    clear tmpd;
    
    % ------ PART IIb: power & cross spectral densities, coherence, powers and peak freq
    % some of the code below needs the peak theta frequency; therefore a minor 
    % variant of WP.nccix is needed (princ channel has to come first, that's all)
    swapi=find(WP.nccix(:,1)==AP.LFPpcInd1 & WP.nccix(:,2)==AP.LFPpcInd1);
    Nccix=WP.nccix;
    Nccix([1 swapi],:)=Nccix([swapi 1],:);
    for j=1:size(Nccix,1)
      ci1=AP.LFPccInd(Nccix(j,1));                  
      ci2=AP.LFPccInd(Nccix(j,2));                
      disp([r(i).segmentType ': (cross-) spectral density (' rawCh(AP.LFPInd(Nccix(j,1))).nm ',' rawCh(AP.LFPInd(Nccix(j,2))).nm ')']);
      % ------ compute power density & scale:
      % - multiply by 2 because we need the one-sided power spectrum
      % - multiply by scaling factor for window, see above
      % - divide by number of points in original (NOT zero-padded) data
      % segment
      % - divide by sampling freq to obtain spectral DENSITY
      % - note that cross spectral densities are complex, power spectral densities aren't
      % - the product of a complex number and its conjugate is the same as its absolute
      % value (=magnitude) squared
      % could possibly include xx% confidence values here - see code in csd.m  
      % and supplementary material of Berg et al Science 2007
      tmpP=(tmpDFT(:,:,Nccix(j,1)).*conj(tmpDFT(:,:,Nccix(j,2))))*scalfac*2/AP.ppSeg*WP.osi*1e-6;      
      r(i).rawPMn{ci1,ci2}=mean(tmpP,2);
      r(i).rawPStd{ci1,ci2}=std(tmpP,1,2);
      % comodulogram: only for autospectra
      if ci1==ci2
        for fbix=1:length(comFStartIx)
          % while downsampling spectrum reshape such that each row is an
          % observation (segment) and each column a variable (freq
          tmpRedP(:,fbix)=mean(tmpP(comFStartIx(fbix):comFStartIx(fbix)+comNFBbin-1,:),1)';
        end        
        [r(i).fComod(:,:,Nccix(j,1)),r(i).fComodP(:,:,Nccix(j,1))]=corrcoef(tmpRedP);
      end
      % for computation of powers, peak frequencies, spectral edge
      % frequencies and gamma centroids
      % 1. eliminate line pickup = interpolate at frequencies n*60+/-1 bin
      % 2. use the absolute, rather than the complex-valued spectral densities for segment-wise computations
      % 3. but: retain complex-valued original
      for h=1:length(ix1)
        % this substitutes 60 Hz peaks and harmonics by mean of neighboring bins
        tmpP(ix1{h},:)=repmat(mean(tmpP(ix2{h},:),1),[3 1]);          
      end
      % absolute values in tmpPa are needed for spectral edge freq here 
      % as well as further down for the computation of power in freq bands
      tmpPa=abs(tmpP);      
      % integral of psd
      tmpIntPsd=cumsum(tmpPa);
      % normalized
      tmpIntPsd=tmpIntPsd./repmat(tmpIntPsd(end,:),size(tmpIntPsd,1),1);
      % compute 95 % SEF
      [tmpNada,tmpEFIx]=min(abs(tmpIntPsd-.95));
      r(i).rawPSEF{ci1,ci2}=F(tmpEFIx);
      r(i).rawPSEFMn{ci1,ci2}=mean(r(i).rawPSEF{ci1,ci2});
      r(i).rawPSEFStd{ci1,ci2}=std(r(i).rawPSEF{ci1,ci2});      
      
      % ----- theta peak freq, running:
      for k=r(i).ni:-1:1
        % multiplying each of the DFTs with the filter's DFT and computing the cross
        % spectral density is the same as multiplying the cross spectral density by the
        % abs squared of that filter
        mnP=abs(tmpP(:,k)).*bftP;
        % smooth spectrum by convoluting with triangular window
        mnPs=conv(mnP(pPeakix),krnl);
        % get rid of krnlhPts border points on either end
        mnPs=mnPs(1+krnlhPts:end-krnlhPts);
        % and since integration requires inthPts on either side of the peaks,
        % look for peaks only in restricted range of mnPs
        tmpr=evdeal(mnPs(1+inthPts:end-inthPts),'idx','allpeaks');
        % offset-corrected (relative to mnPs) positive peak freq
        pp=tmpr.posPeakT{1}+inthPts;
        % if no peak found, entry for that segment will remain nan (as preset in ccRunTemplate)
        if isfinite(pp)
          integr=[];
          for pix=1:length(pp)
            % among the peaks found pick the most prominent one: go a few bins to left
            % and right of the peaks and calculate area
            integr(pix)=sum(mnPs([-inthPts:inthPts]+pp(pix)));
          end
          % pick peak with largest integral..
          [m,ix]=max(integr);
          % ..not forgetting to add time offset: pPeakix AND krnlhPts due to convolution
          % (=extension by krnlhPts on either side) & subsequent cutting off of 2*krnlPts
          % on either side
          r(i).rawPPeakT{ci1,ci2}(k)=F(pp(ix)+krnlhPts-1+pPeakix(1)-1);
          r(i).rawPPeak{ci1,ci2}(k)=tmpr.posPeak{1}(ix);
        end
      end
      % mean & std of running theta peak freq
      OKix=isfinite(r(i).rawPPeak{ci1,ci2});
      r(i).rawPPeakMn{ci1,ci2}=mean(r(i).rawPPeak{ci1,ci2}(OKix));
      r(i).rawPPeakStd{ci1,ci2}=std(r(i).rawPPeak{ci1,ci2}(OKix));
      r(i).rawPPeakTMn{ci1,ci2}=mean(r(i).rawPPeakT{ci1,ci2}(OKix));
      r(i).rawPPeakTStd{ci1,ci2}=std(r(i).rawPPeakT{ci1,ci2}(OKix));

      % ----- peak freq, derived from average:
      % multiplying each of the DFTs with the filter's DFT and computing the cross
      % spectral density is the same as multiplying the cross spectral density by the 
      % abs squared of that filter
      % NOTE: it is important to FIRST average and THEN compute absolute value
      mnP=abs(mean(tmpP,2)).*bftP;
      % smooth spectrum by convoluting with triangular window
      mnPs=conv(mnP(pPeakix),krnl);
      % get rid of krnlhPts border points on either end
      mnPs=mnPs(krnlhPts+1:end-krnlhPts);            
      % and since integration requires inthPts on either side of the peaks,
      % look for peaks only in restricted range of mnPs
      tmpr=evdeal(mnPs(1+inthPts:end-inthPts),'idx','allpeaks');
      % offset-corrected (relative to mnPs) positive peak freq
      pp=tmpr.posPeakT{1}+inthPts;
      % if no peak found, entry for that channel will remain nan (default entry of
      % WP.ccDerTemplate)
      if isfinite(pp) 
        integr=[];
        for pix=1:length(pp)
          % among the peaks found pick the most prominent one: go a few bins to left
          % and right of the peaks and calculate area 
          integr(pix)=sum(mnPs([-inthPts:inthPts]+pp(pix)));
        end
        % pick peak with largest integral..
        [m,ix]=max(integr);
        % ..not forgetting to add time offset: pPeakix AND krnlhPts due to convolution 
        % (=extension by krnlhPts on either side) & subsequent cutting off of 2*krnlPts 
        % on either side 
        r(i).rawPMnPeakT{ci1,ci2}=F(pp(ix)+krnlhPts-1+pPeakix(1)-1);
        r(i).rawPMnPeak{ci1,ci2}=tmpr.posPeak{1}(ix);
      end
      % ----- plot 2:
      % original spec - should be equivalent to abs(mean(tmpP,2)) (except for line pickup)
      subplot(2,1,2), cla, hold on
      ph(1)=plot(F(prngix),abs(r(i).rawPMn{ci1,ci2}(prngix)),'k');            
      % windowed spec
      ph(2)=plot(F(prngix),mnP(prngix),'b');            
      % windowed & smoothed: 
      ph(3)=plot(F(pPeakix),mnPs,'g');            
      % peak derived from averaged spectra
      php=plot(r(i).rawPMnPeakT{ci1,ci2},r(i).rawPMnPeak{ci1,ci2},'ro');
      set(php,'markerfacecolor','r');
      % average of peaks from individual spectra
      php=plot(r(i).rawPPeakTMn{ci1,ci2},r(i).rawPPeakMn{ci1,ci2},'kd');
      set(php,'markerfacecolor','m');
      errorcross([r(i).rawPPeakTMn{ci1,ci2} r(i).rawPPeakMn{ci1,ci2}],...
        [r(i).rawPPeakTStd{ci1,ci2} r(i).rawPPeakStd{ci1,ci2}],...
        'color','m');
      niceyax;
      title(['cross/power spectra: ' r(i).segmentType ', ' rawCh(AP.LFPInd(Nccix(j,1))).nm '-' rawCh(AP.LFPInd(Nccix(j,2))).nm]);
      xlabel('freq (Hz)');
      ylabel('power spectral density (mV^2/Hz)');
      legend(ph,{'raw','windowed','smoothed'});
      drawnow
      % *** while we're marvelling at the peaks & valleys the CPU should be kept busy 
      % computing power and coherence in freq bands ***
      % For the narrow variant of theta use frequency band (+/- 1 Hz) centered
      % at the theta peak in principal channel at the current behavior
      % If there is no peak we have a pathological case (psd with a slope < 0
      % over whole theta range) rendering the narrow theta parameters
      % meaningless. Use center of theta lo freq range in that case.
      if j==1
        if isnan(r(i).rawPMnPeakT{ci1,ci2})
          logstr{end+1}=['no theta peak found - setting theta narrow range to ' num2str(mean(AP.thetaLo)) ' +/-1 Hz'];
          warning(logstr{end});
          thNarrowFIx=F>=mean(AP.thetaLo)-1 & F<=mean(AP.thetaLo)+1;
        else
          thNarrowFIx=F>=r(i).rawPMnPeakT{ci1,ci2}-1 & F<=r(i).rawPMnPeakT{ci1,ci2}+1;
        end
      end
      % ----- power: integral of psd
      % - segment-wise
      r(i).rawDePE{ci1,ci2}=sum(tmpPa(deFIx,:),1)*fbinw;      
      r(i).rawThPE{ci1,ci2}=sum(tmpPa(thFIx,:),1)*fbinw;
      r(i).rawThNarrowPE{ci1,ci2}=sum(tmpPa(thNarrowFIx,:),1)*fbinw;
      r(i).rawBePE{ci1,ci2}=sum(tmpPa(beFIx,:),1)*fbinw;      
      r(i).rawGaPE{ci1,ci2}=sum(tmpPa(gaFIx,:),1)*fbinw;
      r(i).rawGaNarrowPE{ci1,ci2}=sum(tmpPa(gaNarrowFIx,:),1)*fbinw;      
      r(i).rawRiPE{ci1,ci2}=sum(tmpPa(riFIx,:),1)*fbinw;      
      % - means & std
      r(i).rawDePEMn{ci1,ci2}=mean(r(i).rawDePE{ci1,ci2});
      r(i).rawDePEStd{ci1,ci2}=std(r(i).rawDePE{ci1,ci2});          
      r(i).rawThPEMn{ci1,ci2}=mean(r(i).rawThPE{ci1,ci2});
      r(i).rawThPEStd{ci1,ci2}=std(r(i).rawThPE{ci1,ci2});          
      r(i).rawThNarrowPEMn{ci1,ci2}=mean(r(i).rawThNarrowPE{ci1,ci2});
      r(i).rawThNarrowPEStd{ci1,ci2}=std(r(i).rawThNarrowPE{ci1,ci2});          
      r(i).rawBePEMn{ci1,ci2}=mean(r(i).rawBePE{ci1,ci2});
      r(i).rawBePEStd{ci1,ci2}=std(r(i).rawBePE{ci1,ci2});          
      r(i).rawGaPEMn{ci1,ci2}=mean(r(i).rawGaPE{ci1,ci2});
      r(i).rawGaPEStd{ci1,ci2}=std(r(i).rawGaPE{ci1,ci2});          
      r(i).rawGaNarrowPEMn{ci1,ci2}=mean(r(i).rawGaNarrowPE{ci1,ci2});
      r(i).rawGaNarrowPEStd{ci1,ci2}=std(r(i).rawGaNarrowPE{ci1,ci2});          
      r(i).rawRiPEMn{ci1,ci2}=mean(r(i).rawRiPE{ci1,ci2});
      r(i).rawRiPEStd{ci1,ci2}=std(r(i).rawRiPE{ci1,ci2});          
%       % - from average
%       r(i).rawPMnDeP{ci1,ci2}=sum(abs(r(i).rawPMn{ci1,ci2}(deFIx)))*fbinw;
%       r(i).rawPMnThP{ci1,ci2}=sum(abs(r(i).rawPMn{ci1,ci2}(thFIx)))*fbinw;
%       r(i).rawPMnThNarrowP{ci1,ci2}=sum(abs(r(i).rawPMn{ci1,ci2}(thNarrowFIx)))*fbinw;
%       r(i).rawPMnBeP{ci1,ci2}=sum(abs(r(i).rawPMn{ci1,ci2}(beFIx)))*fbinw;
%       r(i).rawPMnGaP{ci1,ci2}=sum(abs(r(i).rawPMn{ci1,ci2}(gaFIx)))*fbinw;
%       r(i).rawPMnRiP{ci1,ci2}=sum(abs(r(i).rawPMn{ci1,ci2}(riFIx)))*fbinw;

      % ----- gamma centroids: frequencies weighted by power
      % - segment-wise
      r(i).rawGaCentroid{ci1,ci2}=sum(tmpPa(gaFIx,:).*repmat(F(gaFIx),1,r(i).ni),1)./sum(tmpPa(gaFIx,:),1);      
      % - means & std
      r(i).rawGaCentroidMn{ci1,ci2}=mean(r(i).rawGaCentroid{ci1,ci2});
      r(i).rawGaCentroidStd{ci1,ci2}=std(r(i).rawGaCentroid{ci1,ci2});
      
      % coherence
      if ci1==ci2
        r(i).rawCohMn{ci1,ci2}=ones([AP.ppSeg/2 1]);
        r(i).rawCohStd{ci1,ci2}=zeros([AP.ppSeg/2 1]);        
        % coherence in freq bands
        r(i).rawCohMnDe{ci1,ci2}=1;
        r(i).rawCohStdDe{ci1,ci2}=0;
        r(i).rawCohMnTh{ci1,ci2}=1;
        r(i).rawCohStdTh{ci1,ci2}=0;
        r(i).rawCohMnThNarrow{ci1,ci2}=1;
        r(i).rawCohStdThNarrow{ci1,ci2}=0;
        r(i).rawCohMnBe{ci1,ci2}=1;
        r(i).rawCohStdBe{ci1,ci2}=0;
        r(i).rawCohMnGa{ci1,ci2}=1;
        r(i).rawCohStdGa{ci1,ci2}=0;
        r(i).rawCohMnGaNarrow{ci1,ci2}=1;
        r(i).rawCohStdGaNarrow{ci1,ci2}=0;
        r(i).rawCohMnRi{ci1,ci2}=1;        
        r(i).rawCohStdRi{ci1,ci2}=0;
      else
        r(i).rawCohMn{ci1,ci2}=abs(r(i).rawPMn{ci1,ci2}).^2./(r(i).rawPMn{ci1,ci1}.*r(i).rawPMn{ci2,ci2});
        r(i).rawCohStd{ci1,ci2}=abs(r(i).rawPStd{ci1,ci2}).^2./(r(i).rawPStd{ci1,ci1}.*r(i).rawPStd{ci2,ci2});
        % coherence in freq bands:
        r(i).rawCohMnDe{ci1,ci2}=mean(r(i).rawCohMn{ci1,ci2}(deFIx));
        r(i).rawCohStdDe{ci1,ci2}=mean(r(i).rawCohStd{ci1,ci2}(deFIx));
        r(i).rawCohMnTh{ci1,ci2}=mean(r(i).rawCohMn{ci1,ci2}(thFIx));
        r(i).rawCohStdTh{ci1,ci2}=mean(r(i).rawCohStd{ci1,ci2}(thFIx));
        r(i).rawCohMnThNarrow{ci1,ci2}=mean(r(i).rawCohMn{ci1,ci2}(thNarrowFIx));
        r(i).rawCohStdThNarrow{ci1,ci2}=mean(r(i).rawCohStd{ci1,ci2}(thNarrowFIx));
        r(i).rawCohMnBe{ci1,ci2}=mean(r(i).rawCohMn{ci1,ci2}(beFIx));
        r(i).rawCohStdBe{ci1,ci2}=mean(r(i).rawCohStd{ci1,ci2}(beFIx));
        r(i).rawCohMnGa{ci1,ci2}=mean(r(i).rawCohMn{ci1,ci2}(gaFIx));
        r(i).rawCohStdGa{ci1,ci2}=mean(r(i).rawCohStd{ci1,ci2}(gaFIx));
        r(i).rawCohMnGaNarrow{ci1,ci2}=mean(r(i).rawCohMn{ci1,ci2}(gaNarrowFIx));
        r(i).rawCohStdGaNarrow{ci1,ci2}=mean(r(i).rawCohStd{ci1,ci2}(gaNarrowFIx));
        r(i).rawCohMnRi{ci1,ci2}=mean(r(i).rawCohMn{ci1,ci2}(riFIx));
        r(i).rawCohStdRi{ci1,ci2}=mean(r(i).rawCohStd{ci1,ci2}(riFIx));        
      end
      % pause;
    end % for:size(Nccix,1)
  end % if: not isempty segmentType
end
close(fh_spDet);      
% bring up the summary results figure
drawnow;




% appendix: coherence for dummies (like the author of this m-file)
% Sth. like a segment-wise ('running') coherence cannot be computed (if we stick 
% to the precise definition of coherence). Why? If we denote data segments from two 
% channels as x and y, coherence is 
%         abs(cross spectrum(x,y)).^2 / (power spec(x) .* power spec(y))
% so, denoting X and Y as the discrete fourier transforms of x and y
%         abs(X .* conj(Y)).^2 / (X .* conj(X) * Y .* conj(Y))
% which will ALWAYS yield 1 if x and y are single segments. Put simply, it is 
%         (x1*y1)^2/((x1*x1) * (y1*y1)) = 1
% 
% If X and Y are averaged discrete fourier transforms the numerator will always be leq 
% the denominator, resulting in coherence values between 0 and 1.
%
% x1=5;
% x2=3;
% y1=4;
% y2=7;
% 
% (x1*y1 + x2*y2).^2/((x1.^2 + x2.^2) .* (y1.^2 + y2.^2))