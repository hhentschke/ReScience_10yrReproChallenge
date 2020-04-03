function rmouse_cspecp_envelope(rawCh,strmType,STshort)
% This is a stripped-down version of rmouse_cspecP dealing with 
% envelope streams.

global DS AP WP r logstr


% ------ PART I: preparations
% check whether AP.ppSeg is really a power of 2
if log2(AP.ppSeg)~=round(log2(AP.ppSeg)), 
  error('AP.ppSeg must be a power of 2'); 
end
% window
if strcmpi(AP.dftWin,'rect'), window=ones(AP.ppSeg,1); 
else eval(['window=' AP.dftWin '(AP.ppSeg);']);
end

% scale the window: sum of squared elements of window
% divided by (scalfac*AP.ppSeg) must be==1
scalfac=AP.ppSeg/sum(window.^2);

% decide to what frequency range to limit the results
switch strmType
  case 'gammaEnv'
    % for any kind of stream a spectrum of its envelope does not make sense
    % for frequencies above its lower corner freq. However, for plots
    % intended to show exactly this, go up to upper corner freq.
    fLim=[0 AP.gamma(2)];
  case 'gammaNarrowEnv'
    fLim=[0 AP.gammaNarrow(2)];
  otherwise
    error('illegal stream type in spectral analysis of streams');
end

% frequency: only in 1st element of r (is same for all segmentTypes and streams)
F=makecol(1e6/WP.osi/AP.ppSeg*[0:AP.ppSeg/2-1]);
fbinw=diff(F(1:2));
% now restict
fLimIx=find(F>=fLim(1) & F<=fLim(2));
F=F(fLimIx);
eval(['r(1).' STshort 'F=F;']);

% half width, in bins, of triangular kernel for smoothing of spectra 
% (more exactly, kernel will have length 2*krnlhPts+1)
krnlhPts=ceil(AP.specKrnlHWid/fbinw);
krnl=triang(2*krnlhPts+1);
% dont forget to normalize
krnl=krnl/sum(krnl);
% half width of interval for computation of peak integral in bins
inthPts=round(AP.peakHW/fbinw);
% indices into F (freq bands) (limit to delta and theta bands for envelope
% streams)
deFIx=F>=AP.delta(1) & F<=AP.delta(2);
thFIx=F>=AP.theta(1) & F<=AP.theta(2);

prngix=find(F>=min(AP.peakFRng(1),AP.FIRcf(1)) & F<=max(AP.peakFRng(2),AP.FIRcf(2)));      

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
bft=fft(fir1(AP.FIRfo,AP.FIRcf*1e-6*WP.osi*2),AP.ppSeg);
% keep one-sided DFT only and make a column vector
bft=(bft(1:AP.ppSeg/2))';
% restrict freq range
bft=bft(fLimIx);
% 'unscaled' power is also needed:
bftP=abs(bft).^2;

% --- figure 
tmpftag='SpecPeakDet_envelope';
fh_spDet=findobj('tag',tmpftag);
if isempty(fh_spDet), fh_spDet=figure;
else  figure(fh_spDet);
end
tmpScrSz=get(0,'Screensize');
tmpScrSz([3 4])=round([tmpScrSz(3)*.35 tmpScrSz(4)*.8]);  
set(fh_spDet,'position',tmpScrSz,'tag',tmpftag,'name','Peak detection in envelope stream power spectra' ,...
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
    ccRunTemplate=WP.ccDerTemplate;
    % ..initialized with vector of nans
    ccRunTemplate(AP.trixie)={tmpvon};
    % assign preallocated var templates:
    % cross and power spectral density (whole freq range): averages
    strm__PMn=WP.psTemplate;
    strm__PStd=WP.psTemplate;
    % power in freq bands: running & averages 
    strm__DePE=ccRunTemplate;
    strm__DePEMn=WP.ccDerTemplate;
    strm__DePEStd=WP.ccDerTemplate;
    strm__ThPE=ccRunTemplate;
    strm__ThPEMn=WP.ccDerTemplate;
    strm__ThPEStd=WP.ccDerTemplate;
    
    strm__ThNarrowPE=ccRunTemplate;
    strm__ThNarrowPEMn=WP.ccDerTemplate;
    strm__ThNarrowPEStd=WP.ccDerTemplate;

    % coherence (whole freq range): averages
    strm__CohMn=WP.psTemplate;
    strm__CohStd=WP.psTemplate;
    % average coherence in frequency bands: derived from averages (see bottom of mfile)
    strm__CohMnDe=WP.ccDerTemplate;
    strm__CohStdDe=WP.ccDerTemplate;    
    strm__CohMnTh=WP.ccDerTemplate;
    strm__CohStdTh=WP.ccDerTemplate;    
    strm__CohMnThNarrow=WP.ccDerTemplate;
    strm__CohStdThNarrow=WP.ccDerTemplate;    
    % theta peak magnitude and freq: derived from average
    strm__PMnPeak=WP.ccDerTemplate;
    strm__PMnPeakT=WP.ccDerTemplate;    
    % theta peak magnitude and freq: running & averages
    strm__PPeak=ccRunTemplate;
    strm__PPeakT=ccRunTemplate;
    strm__PPeakMn=WP.ccDerTemplate;
    strm__PPeakStd=WP.ccDerTemplate;
    strm__PPeakTMn=WP.ccDerTemplate;    
    strm__PPeakTStd=WP.ccDerTemplate;    
    % complex intermediate results var holding DFT from all LFP channels to be analyzed
    tmpDFT=repmat(nan+sqrt(-1),[length(fLimIx) r(i).ni AP.nLFPCh ]);
    % disp(['intermediate results var takes up ' num2str(prod(size(tmpDFT))*8/2^20,'%4.1f') ' Mb']);            
    % ------ PART IIa: DFT of each LFP channel
    for chInd=1:AP.nLFPCh 
      disp([r(i).segmentType ': DFT ' strmType ' ch ' rawCh(AP.LFPInd(chInd)).nm ' ..']);
      % load data
      eval(['tmpd=strmread([AP.strmDir ''\'' rawCh(AP.LFPInd(chInd)).' strmType 'Fn],''nPts'',r(i).iPts(end,2),''verbose'',0);']);
      for k=r(i).ni:-1:1
        tmpIdx=r(i).iPts(k,1):r(i).iPts(k,2);
        % data segments are detrended individually
        tmpfft=fft(detrend(tmpd(tmpIdx),'constant').*window);
        % keep one-sided
        tmpDFT(:,k,chInd)=tmpfft(fLimIx);
      end
    end
    clear tmpd;
    
    % ------ PART IIb: power spectral densities, coherence and peak freq
    % the part computing average coherence at theta needs the peak theta frequency;
    % therefore a minor variant of WP.nccix is needed (princ channel has to come first,
    % that's all)
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
      % - divide by number of points in original (NOT zero-padded) data segment
      % - divide by sampling freq to obtain spectral DENSITY
      % - note that cross spectral densities are complex, power spectral densities aren't
      % - the product of a complex number and its conjugate is the same as its absolute
      % value (=magnitude) squared
      % could possibly include xx% confidence values here - see code in csd.m      
      tmpP=(tmpDFT(:,:,Nccix(j,1)).*conj(tmpDFT(:,:,Nccix(j,2))))*scalfac*2/AP.ppSeg*WP.osi*1e-6;      
      strm__PMn{ci1,ci2}=mean(tmpP,2);
      strm__PStd{ci1,ci2}=std(tmpP,1,2);

      % for computation of powers and peak frequencies 
      % 2. use the absolute, rather than the complex-valued spectral densities for segment-wise computations
      % 3. but: retain complex-valued original
      % (tmpPa will be needed further down for the computation of power in freq
      % bands)
      tmpPa=abs(tmpP);      
      
      % ----- peak freq, running:
      for k=r(i).ni:-1:1
        % multiplying each of the DFTs with the filter's DFT and computing the cross
        % spectral density is the same as multiplying the cross spectral density by the
        % abs squared of that filter
        mnP=abs(tmpP(:,k)).*bftP;
        % smooth spectrum by convoluting with triangular window
        mnPs=conv(mnP(pPeakix),krnl);
        % get rid of krnlhPts border points on either end
        mnPs=mnPs(krnlhPts+1:end-krnlhPts);
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
          strm__PPeakT{ci1,ci2}(k)=F(pp(ix)+krnlhPts-1+pPeakix(1)-1);
          strm__PPeak{ci1,ci2}(k)=tmpr.posPeak{1}(ix);
        end
      end
      % mean & std 
      OKix=isfinite(strm__PPeak{ci1,ci2});
      strm__PPeakMn{ci1,ci2}=mean(strm__PPeak{ci1,ci2}(OKix));
      strm__PPeakStd{ci1,ci2}=std(strm__PPeak{ci1,ci2}(OKix));
      strm__PPeakTMn{ci1,ci2}=mean(strm__PPeakT{ci1,ci2}(OKix));
      strm__PPeakTStd{ci1,ci2}=std(strm__PPeakT{ci1,ci2}(OKix));

      % ----- peak freq, derived from average:
      % multiplying each of the DFTs with the filter's DFT and computing the cross
      % spectral density is the same as multiplying the cross spectral density by the 
      % abs squared of that filter
      % NOTE: it is important to FIRST average and THEN compute absolute value
      mnP=abs(mean(tmpP,2)).*bftP;
      % smooth spectrum by convoluting with triangular window
      mnPs=conv(mnP(pPeakix),krnl);
      % get rid of krnlhPts border points on either end
      mnPs=mnPs(1+krnlhPts:end-krnlhPts);
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
        strm__PMnPeakT{ci1,ci2}=F(pp(ix)+krnlhPts-1+pPeakix(1)-1);
        strm__PMnPeak{ci1,ci2}=tmpr.posPeak{1}(ix);
      end
      % ----- plot 2:
      % original spec - should be equivalent to abs(mean(tmpP,2)) (except for line pickup)
      subplot(2,1,2), cla, hold on
      ph(1)=plot(F(prngix),abs(strm__PMn{ci1,ci2}(prngix)),'k');            
      % windowed spec
      ph(2)=plot(F(prngix),mnP(prngix),'b');            
      % windowed & smoothed: 
      ph(3)=plot(F(pPeakix),mnPs,'g');            
      % peak derived from averaged spectra
      php=plot(strm__PMnPeakT{ci1,ci2},strm__PMnPeak{ci1,ci2},'ro');
      set(php,'markerfacecolor','r');
      % average of peaks from individual spectra
      php=plot(strm__PPeakTMn{ci1,ci2},strm__PPeakMn{ci1,ci2},'kd');
      set(php,'markerfacecolor','m');
      errorcross([strm__PPeakTMn{ci1,ci2} strm__PPeakMn{ci1,ci2}],...
        [strm__PPeakTStd{ci1,ci2} strm__PPeakStd{ci1,ci2}],...
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
      % at the peak in principal channel at the current behavior
      % If there is no peak we have a pathological case (psd with a slope < 0
      % over whole theta range) rendering the narrow theta parameters
      % meaningless. Use center of theta lo freq range in that case.
      if j==1
        if isnan(strm__PMnPeakT{ci1,ci2})
          logstr{end+1}=['no theta peak found - setting theta narrow range to ' num2str(mean(AP.thetaLo)) ' +/-1 Hz'];
          thNarrowFIx=F>=mean(AP.thetaLo)-1 & F<=mean(AP.thetaLo)+1;
        else
          thNarrowFIx=F>=strm__PMnPeakT{ci1,ci2}-1 & F<=strm__PMnPeakT{ci1,ci2}+1;
        end
      end
      % ----- power: integral of psd
      % - segment-wise
      strm__DePE{ci1,ci2}=sum(tmpPa(deFIx,:),1)*fbinw;
      strm__ThPE{ci1,ci2}=sum(tmpPa(thFIx,:),1)*fbinw;
      strm__ThNarrowPE{ci1,ci2}=sum(tmpPa(thNarrowFIx,:),1)*fbinw;      
      % - means & std
      strm__DePEMn{ci1,ci2}=mean(strm__DePE{ci1,ci2});
      strm__DePEStd{ci1,ci2}=std(strm__DePE{ci1,ci2});
      strm__ThPEMn{ci1,ci2}=mean(strm__ThPE{ci1,ci2});
      strm__ThPEStd{ci1,ci2}=std(strm__ThPE{ci1,ci2});
      strm__ThNarrowPEMn{ci1,ci2}=mean(strm__ThNarrowPE{ci1,ci2});
      strm__ThNarrowPEStd{ci1,ci2}=std(strm__ThNarrowPE{ci1,ci2});
      
      if ci1==ci2
        strm__CohMn{ci1,ci2}=ones([length(fLimIx) 1]);
        strm__CohStd{ci1,ci2}=zeros([length(fLimIx) 1]);
        % coherence in freq bands
        strm__CohMnDe{ci1,ci2}=1;
        strm__CohStdDe{ci1,ci2}=0;
        strm__CohMnTh{ci1,ci2}=1;
        strm__CohStdTh{ci1,ci2}=0;
        strm__CohMnThNarrow{ci1,ci2}=1;
        strm__CohStdThNarrow{ci1,ci2}=0;
      else
        strm__CohMn{ci1,ci2}=abs(strm__PMn{ci1,ci2}).^2./(strm__PMn{ci1,ci1}.*strm__PMn{ci2,ci2});
        strm__CohStd{ci1,ci2}=abs(strm__PStd{ci1,ci2}).^2./(strm__PStd{ci1,ci1}.*strm__PStd{ci2,ci2});
        % coherence in freq bands:
        strm__CohMnDe{ci1,ci2}=mean(strm__CohMn{ci1,ci2}(deFIx));
        strm__CohStdDe{ci1,ci2}=mean(strm__CohStd{ci1,ci2}(deFIx));
        strm__CohMnTh{ci1,ci2}=mean(strm__CohMn{ci1,ci2}(thFIx));
        strm__CohStdTh{ci1,ci2}=mean(strm__CohStd{ci1,ci2}(thFIx));
        strm__CohMnThNarrow{ci1,ci2}=mean(strm__CohMn{ci1,ci2}(thNarrowFIx));
        strm__CohStdThNarrow{ci1,ci2}=mean(strm__CohStd{ci1,ci2}(thNarrowFIx));
      end
      % pause;
    end % for:size(Nccix,1)

    % now, at the end, give computed results proper name and make fields of r
    tmpSVar=whos('strm__*');
    for tmpVix=1:length(tmpSVar)
      eval(['r(i).' STshort tmpSVar(tmpVix).name(7:end) '=' tmpSVar(tmpVix).name '; clear ' tmpSVar(tmpVix).name]);  
    end
    
  end % if: not isempty segmentType
end
close(fh_spDet);      


