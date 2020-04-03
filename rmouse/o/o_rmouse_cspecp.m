% computes power and cross power spectral density
% to do:
% - update computation of peak memory demand
% ------ PART I: preparations
% check whether AP.ppSeg is really a power of 2
if log2(AP.ppSeg)~=round(log2(AP.ppSeg)), error('AP.ppSeg must be a power of 2'); end
% window
if strcmpi(AP.dftWin,'rect'), window=ones(AP.ppSeg,1); 
else eval(['window=' AP.dftWin '(AP.ppSeg);']);
end
% scale the window: sum of squared elements of window 
% divided by (scalfac*AP.ppSeg) must be==1
scalfac=AP.ppSeg/sum(window.^2);
% frequency: only in 1st element of r (is same for all segmentTypes)
F=makecol(1e6/osi/AP.ppSeg*[0:AP.ppSeg/2-1]);
r(1).F=F;
fbinw=diff(F(1:2));
% half width, in bins, of triangular kernel for smoothing of spectra 
% (more exactly, kernel will have length 2*krnlhPts+1)
krnlhPts=ceil(AP.specKrnlHWid/fbinw);
krnl=triang(2*krnlhPts+1);
% dont forget to normalize
krnl=krnl/sum(krnl);
% half width of interval for computation of peak integral in bins
inthPts=round(AP.peakHW/fbinw);
% indices into F (freq bands)
deix=F>=AP.delta(1) & F<=AP.delta(2);
thix=F>=AP.theta(1) & F<=AP.theta(2);
gaix=F>=AP.gamma(1) & F<=AP.gamma(2);
riix=F>=AP.ripple(1) & F<=AP.ripple(2);      
prngix=F>=AP.peakFRng(1) & F<=AP.peakFRng(2);      
% indices to range for peak detection:
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
% peak detection: DFT of the FIR filter with which to multiply averaged DFT of data.
% length=order of filter must be < AP.ppSeg (otherwise its DFT will be computed
% from truncated version) and should in general be fixed
bft=fft(fir1(AP.FIRfo,AP.peakFRng*1e-6*osi*2),AP.ppSeg);
% keep one-sided DFT only and make a column vector
bft=(bft(1:AP.ppSeg/2))';
% 'unscaled' power is also needed:
bftP=abs(bft).^2;
% preparations for cutting out line hum
humF=60:60:F(end);
for h=1:length(humF)
  % indices to immediately adjacent bins
  [nix,ix]=min(abs(F-humF(h)));
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

% ------- PART II: the works
for i=1:length(r)
  if ~isempty(r(i).iPts)
    disp([r(i).segmentType '..']);
    % 'vector of nans'
    tmpvon=repmat(nan,1,r(i).ni);
    % template for all running pars..
    ccRunTemplate=ccDerTemplate;
    % ..initialized with vector of nans
    ccRunTemplate(AP.trixie)={tmpvon};
    % assign preallocated var templates:
    % cross and power spectral density (whole freq range): averages
    r(i).rawPMn=psTemplate;
    r(i).rawPStd=psTemplate;
    % power in freq bands: running & averages 
    r(i).rawDePE=ccRunTemplate;
    r(i).rawDePEMn=ccDerTemplate;
    r(i).rawDePEStd=ccDerTemplate;
    r(i).rawThPE=ccRunTemplate;
    r(i).rawThPEMn=ccDerTemplate;
    r(i).rawThPEStd=ccDerTemplate;
    r(i).rawGaPE=ccRunTemplate;
    r(i).rawGaPEMn=ccDerTemplate;
    r(i).rawGaPEStd=ccDerTemplate;
    r(i).rawRiPE=ccRunTemplate;
    r(i).rawRiPEMn=ccDerTemplate;
    r(i).rawRiPEStd=ccDerTemplate;
    % coherence (whole freq range): averages
    r(i).rawCohMn=psTemplate;
    r(i).rawCohStd=psTemplate;    
    % average coherence in frequency bands: derived from averages (see bottom of mfile)
    r(i).rawCohMnDe=ccDerTemplate;
    r(i).rawCohStdDe=ccDerTemplate;    
    r(i).rawCohMnTh=ccDerTemplate;
    r(i).rawCohStdTh=ccDerTemplate;    
    r(i).rawCohMnGa=ccDerTemplate;    
    r(i).rawCohStdGa=ccDerTemplate;    
    r(i).rawCohMnRi=ccDerTemplate;
    r(i).rawCohStdRi=ccDerTemplate;    
    % theta peak magnitude and freq: derived from average
    r(i).rawPMnPeak=ccDerTemplate;
    r(i).rawPMnPeakT=ccDerTemplate;    
    % theta peak magnitude and freq: running & averages
    r(i).rawPPeak=ccRunTemplate;
    r(i).rawPPeakT=ccRunTemplate;
    r(i).rawPPeakMn=ccDerTemplate;
    r(i).rawPPeakStd=ccDerTemplate;
    r(i).rawPPeakTMn=ccDerTemplate;    
    r(i).rawPPeakTStd=ccDerTemplate;    
    % complex intermediate results var holding DFT from all LFP channels 2b analyzed
    nix=repmat(nan+sqrt(-1),[AP.ppSeg/2 r(i).ni nLFPCh]);
    % disp(['intermediate results var takes up ' num2str(prod(size(nix))*8/2^20,'%4.1f') ' Mb']);            
    % ------ PART IIa: DFT of each LFP channel
    for chInd=1:nLFPCh
      disp([r(i).segmentType ': DFT raw ch ' rawCh(AP.LFPInd(chInd)).nm ' ..']);
      % load data
      if exist([dpath '\' abfFn '.mat'],'file')
        tmpd=matDload([dpath '\' abfFn '.mat'],'start',0,'stop',discrete2cont(r(i).iPts(end,2),osi*1e-6,'intv',1),...
          'channels',{rawCh(AP.LFPInd(chInd)).nm},'verbose',verbose);
      else
        tmpd=abfload([dpath '\' abfFn '.abf'],'start',0,'stop',discrete2cont(r(i).iPts(end,2),osi*1e-6,'intv',1),...
          'channels',{rawCh(AP.LFPInd(chInd)).nm},'verbose',verbose);
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
          error('abfload or matDload handles time information poorly');
        end
      end
      for k=r(i).ni:-1:1
        tmpIdx=r(i).iPts(k,1):r(i).iPts(k,2);
        % data segments are detrended individually
        tmpfft=fft(detrend(tmpd(tmpIdx),'constant').*window);
        % keep one-sided
        nix(:,k,chInd)=tmpfft(1:AP.ppSeg/2);
      end
    end
    clear tmpd;
    
    % ------ PART IIb: power & cross spectral densities, coherence, powers and peak freq
    % the part computing average coherence at theta needs the peak theta frequency;
    % therefore a minor variant of nccix is needed (princ channel has to come first,
    % that's all)
    swapi=find(nccix(:,1)==AP.LFPpcInd1 & nccix(:,2)==AP.LFPpcInd1);
    Nccix=nccix;
    Nccix([1 swapi],:)=Nccix([swapi 1],:);
    for j=1:size(Nccix,1)
      ci1=AP.LFPccInd(Nccix(j,1));                  
      ci2=AP.LFPccInd(Nccix(j,2));                
      disp([r(i).segmentType ': (cross-) spectral density (' rawCh(AP.LFPInd(Nccix(j,1))).nm ',' rawCh(AP.LFPInd(Nccix(j,2))).nm ')']);
      % ------ compute power density & scale:
      % - multiply by 2 because we need the one-sided power spectrum
      % - multiply by scaling factor for window, see above
      % - divide by number of points in original (NOT zero-padded) data segment
      % - divide by sampling freq to obtain spectral DENSITY
      % - note that cross spectral densities are complex, power spectral densities aren't
      % - the product of a complex number and its conjugate is the same as its absolute
      % value (=magnitude) squared
      tmpP=(nix(:,:,Nccix(j,1)).*conj(nix(:,:,Nccix(j,2))))*scalfac*2/AP.ppSeg*osi*1e-6;      
      r(i).rawPMn{ci1,ci2}=mean(tmpP,2);
      r(i).rawPStd{ci1,ci2}=std(tmpP,1,2);
      % could possibly include xx% confidence values here - see code in csd.m
      % for computation of powers and peak frequencies 
      % 1. eliminate line pickup = interpolate at frequencies n*60+/-1 bin
      % 2. use the absolute, rather than the complex-valued spectral densities for segment-wise computations
      % 3. but: retain complex-valued original
      for h=1:length(ix1)
        % this substitutes 60 Hz peaks and harmonics by mean of neighboring bins
        tmpP(ix1{h},:)=repmat(mean(tmpP(ix2{h},:),1),[3 1]);          
      end
      tmpPa=abs(tmpP);      
      % ----- power: integral of psd
      % - segment-wise
      r(i).rawDePE{ci1,ci2}=sum(tmpPa(deix,:),1)*fbinw;      
      r(i).rawThPE{ci1,ci2}=sum(tmpPa(thix,:),1)*fbinw;
      r(i).rawGaPE{ci1,ci2}=sum(tmpPa(gaix,:),1)*fbinw;
      r(i).rawRiPE{ci1,ci2}=sum(tmpPa(riix,:),1)*fbinw;      
      % - means & std
      r(i).rawDePEMn{ci1,ci2}=mean(r(i).rawDePE{ci1,ci2});
      r(i).rawDePEStd{ci1,ci2}=std(r(i).rawDePE{ci1,ci2});          
      r(i).rawThPEMn{ci1,ci2}=mean(r(i).rawThPE{ci1,ci2});
      r(i).rawThPEStd{ci1,ci2}=std(r(i).rawThPE{ci1,ci2});          
      r(i).rawGaPEMn{ci1,ci2}=mean(r(i).rawGaPE{ci1,ci2});
      r(i).rawGaPEStd{ci1,ci2}=std(r(i).rawGaPE{ci1,ci2});          
      r(i).rawRiPEMn{ci1,ci2}=mean(r(i).rawRiPE{ci1,ci2});
      r(i).rawRiPEStd{ci1,ci2}=std(r(i).rawRiPE{ci1,ci2});          
      
      % ----- peak freq, running:
      for k=r(i).ni:-1:1
        % multiplying each of the DFTs with the filter's DFT and computing the cross
        % spectral density is the same as multiplying the cross spectral density by the
        % abs squared of that filter
        mnP=abs(tmpP(:,k)).*bftP;
        % smooth spectrum by convoluting with triangular window
        mnPs=conv(mnP(pPeakix),krnl);
        % get rid of 2*krnlhPts border points on either end
        mnPs=mnPs(2*krnlhPts:end-2*krnlhPts+1);
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
      % mean & std 
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
      % get rid of 2*krnlhPts border points on either end
      mnPs=mnPs(2*krnlhPts:end-2*krnlhPts+1);            
      % and since integration requires inthPts on either side of the peaks,
      % look for peaks only in restricted range of mnPs
      tmpr=evdeal(mnPs(1+inthPts:end-inthPts),'idx','allpeaks');
      % offset-corrected (relative to mnPs) positive peak freq
      pp=tmpr.posPeakT{1}+inthPts;
      % if no peak found, entry for that channel will remain nan (default entry of
      % ccDerTemplate)
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
      [tmpidx1,tmpidx2]=embedtrc(length(find(prngix)),find(F(prngix)==r(i).rawPMnPeakT{ci1,ci2}),length(mnPs),pp(ix));
      ph(3)=plot(F(prngix),mnPs(tmpidx2(1):tmpidx2(2)),'g');            
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
      % while we're marvelling at the peaks & valleys the CPU should be kept busy computing coherences.
      % In the case of gamma and ripples average coherence over the frequency bands given in AP, but 
      % for theta use a narrow frequency band (+/- AP.peakHW Hz) centered at the peak in principal channel 
      %  at the current behavior
      if j==1
        narrThix=F>=r(i).rawPMnPeakT{ci1,ci2}-AP.peakHW & F<=r(i).rawPMnPeakT{ci1,ci2}+AP.peakHW;
      end
      if ci1==ci2
        r(i).rawCohMn{ci1,ci2}=ones([AP.ppSeg/2 1]);
        r(i).rawCohStd{ci1,ci2}=zeros([AP.ppSeg/2 1]);        
        % coherence in freq bands
        % # include std
        r(i).rawCohMnDe{ci1,ci2}=1;
        r(i).rawCohStdDe{ci1,ci2}=0;
        r(i).rawCohMnTh{ci1,ci2}=1;
        r(i).rawCohStdTh{ci1,ci2}=0;
        r(i).rawCohMnGa{ci1,ci2}=1;
        r(i).rawCohStdGa{ci1,ci2}=0;
        r(i).rawCohMnRi{ci1,ci2}=1;        
        r(i).rawCohStdRi{ci1,ci2}=0;
      else
        r(i).rawCohMn{ci1,ci2}=abs(r(i).rawPMn{ci1,ci2}).^2./(r(i).rawPMn{ci1,ci1}.*r(i).rawPMn{ci2,ci2});
        r(i).rawCohStd{ci1,ci2}=abs(r(i).rawPStd{ci1,ci2}).^2./(r(i).rawPStd{ci1,ci1}.*r(i).rawPStd{ci2,ci2});
        % coherence in freq bands:
        r(i).rawCohMnDe{ci1,ci2}=mean(r(i).rawCohMn{ci1,ci2}(deix));
        r(i).rawCohStdDe{ci1,ci2}=mean(r(i).rawCohStd{ci1,ci2}(deix));
        r(i).rawCohMnTh{ci1,ci2}=mean(r(i).rawCohMn{ci1,ci2}(narrThix));
        r(i).rawCohStdTh{ci1,ci2}=mean(r(i).rawCohStd{ci1,ci2}(narrThix));
        r(i).rawCohMnGa{ci1,ci2}=mean(r(i).rawCohMn{ci1,ci2}(gaix));
        r(i).rawCohStdGa{ci1,ci2}=mean(r(i).rawCohStd{ci1,ci2}(gaix));
        r(i).rawCohMnRi{ci1,ci2}=mean(r(i).rawCohMn{ci1,ci2}(riix));
        r(i).rawCohStdRi{ci1,ci2}=mean(r(i).rawCohStd{ci1,ci2}(riix));        
      end
      % pause;
    end % for:size(Nccix,1)
  end % if: not isempty segmentType
end
clear tmp* swapi Nccix window nix scalfac P F m ix* fbin* narrThix deix thix gaix riix peakix pPeakix mnP* krnl* bft* prngix ph* fnix integr pp ci* inthPts OKix ccRunTe*
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