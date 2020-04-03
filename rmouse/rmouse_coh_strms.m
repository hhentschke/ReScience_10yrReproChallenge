function rmouse_coh_strms(cohType,rawCh)
% ** function rmouse_coh_strms(cohType,rawCh)
% Computes coherence between streams (currently, raw data and gamma
% envelope). All combinations of channels.

global DS AP WP r logstr
% ------ PART I: preparations
switch cohType
  case 'rawgae'
    strmType={'raw','gammaEnv'};
    % for any kind of stream a spectrum of its envelope does not make sense
    % for frequencies above its lower corner freq
    fLim=[0 AP.gamma(1)];
  otherwise
    error('illegal cohType in between-streams coherence computation');
end

% --- check whether AP.ppSeg is really a power of 2
if log2(AP.ppSeg)~=round(log2(AP.ppSeg))
  error('AP.ppSeg must be a power of 2');
end

% --- window
if strcmpi(AP.dftWin,'rect'), window=ones(AP.ppSeg,1); 
else eval(['window=' AP.dftWin '(AP.ppSeg);']);
end
% scale: sum of squared elements of window  divided by (scalfac*AP.ppSeg) 
% must be==1
scalfac=AP.ppSeg/sum(window.^2);

% frequency: only in 1st element of r (is same for all segmentTypes and streams)
F=makecol(1e6/WP.osi/AP.ppSeg*[0:AP.ppSeg/2-1]);
fbinw=diff(F(1:2));
% now restict
fLimIx=find(F>=fLim(1) & F<=fLim(2));
F=F(fLimIx);
eval(['r(1).' cohType 'F=F;']);

% --- kernel for smoothing of spectra
% half width, in bins, of triangular kernel (more exactly, kernel will 
% have length 2*krnlhPts+1)
krnlhPts=ceil(AP.specKrnlHWid/fbinw);
krnl=triang(2*krnlhPts+1);
% dont forget to normalize
krnl=krnl/sum(krnl);

% --- half width of interval for computation of peak integral in bins
inthPts=round(AP.peakHW/fbinw);

% --- indices into F (freq bands)
thFIx=F>=AP.theta(1) & F<=AP.theta(2);

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
  logstr{end+1}=['interval for peak detection in coherence spectra needed to be adjusted to [' num2str(makerow(F(peakix))) '] Hz'];
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
% restrict freq range
bft=bft(fLimIx);
% 'unscaled' power is also needed:
bftP=abs(bft).^2;


% ------------ PART II: the works
for i=1:length(r)
  if ~isempty(r(i).iPts)
    disp([r(i).segmentType '..']);
    % preallocate results variables  
    % *** note: the results variables strm_Coh* hold data for ALL LFP
    % channels, that is, they contain NaNs for non-analyzed channels.
    % Other variables like tmpDFT* hold data only for analyzed channels.
    % i. 3D array, each column the coherence derived from streams of
    % channel pair
    strm_Coh=repmat(nan,[length(fLimIx) AP.nAllLFPCh AP.nAllLFPCh]);
    % ii. amplitudes and frequency of peak of coherence in theta range
    strm_CohPeak=strm_Coh(1,:,:);
    strm_CohPeakF=strm_Coh(1,:,:);
    % iii. integral of coherence in theta range
    strm_CohTh=strm_Coh(1,:,:);
    % iv. complex intermediate results vars holding DFT of the two stream 
    % types from all LFP channels to be analyzed
    tmpDFT1=repmat(nan+sqrt(-1),[length(fLimIx) r(i).ni AP.nLFPCh]);
    tmpDFT2=tmpDFT1;

    % ------ PART IIa: DFT of each LFP channel
    for chInd=1:AP.nLFPCh
      disp(['DFT ch ' rawCh(AP.LFPInd(chInd)).nm ' ..']);
      % --- load data stream type 1
      if exist([DS.dpath '\' DS.abfFn '.mat'],'file')
        tmpd1=matDload([DS.dpath '\' DS.abfFn '.mat'],'start',0,'stop',discrete2cont(r(i).iPts(end,2),WP.osi*1e-6,'intv',1),...
          'channels',{rawCh(AP.LFPInd(chInd)).nm},'verbose',0);
      elseif exist([DS.dpath '\' DS.abfFn '.raw'],'file')
        rawFInfo=rawinfo('filename',[DS.dpath '\' DS.abfFn '.raw'],'print','no');
        tmpd1=rawload([DS.dpath '\' DS.abfFn '.raw'],{rawCh(AP.LFPInd(chInd)).nm},...
          [0 discrete2cont(r(i).iPts(end,2),WP.osi*1e-3,'intv',1)],rawFInfo);
        % put into array and convert to mV
        tmpd1=cat(2,tmpd1{:})/1000;
        tmp=rawfi2genfi(rawFInfo,'channels',AP.rawChAnNm);
        si=tmp.si;
      else
        tmpd1=abfload([DS.dpath '\' DS.abfFn '.abf'],'start',0,'stop',discrete2cont(r(i).iPts(end,2),WP.osi*1e-6,'intv',1),...
          'channels',{rawCh(AP.LFPInd(chInd)).nm},'verbose',0);
      end
      if DS.rawSignalInverted
        tmpd1=-1*tmpd1;
      end
      % Due to the fact that abfload computes discrete time from continuous time a little
      % differently than is done in rmouse_xx, in some cases a single data point is
      % amiss. This is bad, and calls for a revision of abfload. In the
      % meantime, here's the workaround:
      tmppd=r(i).iPts(end,2)-length(tmpd1);
      if tmppd,
        if tmppd==1
          disp('missed a point');
          tmpd1(end+1)=tmpd1(end);
        elseif tmppd== -1
          disp('one point too much');
        else
          error('abfload, matDload or rawload handles time information poorly');
        end
      end

      % --- load data stream type 2
      eval(['tmpd2=strmread([AP.strmDir ''\'' rawCh(AP.LFPInd(chInd)).' ...
        strmType{2} 'Fn],''nPts'',r(i).iPts(end,2),''verbose'',0);']);

      for k=r(i).ni:-1:1
        tmpIdx=r(i).iPts(k,1):r(i).iPts(k,2);
        % fft stream 1
        tmpfft=fft(detrend(tmpd1(tmpIdx),'constant').*window);
        % keep only freq values of interest, one-sided
        tmpDFT1(:,k,chInd)=tmpfft(fLimIx);
        % fft stream 2
        tmpfft=fft(detrend(tmpd2(tmpIdx),'constant').*window);
        % keep only freq values of interest, one-sided
        tmpDFT2(:,k,chInd)=tmpfft(fLimIx);
      end
    end
    clear tmpd*

    % ------ PART IIb: power spectral densities, coherence and peak freq
    % ------ how to compute power densities & scale:
    % - multiply by 2 because we need the one-sided power spectrum
    % - multiply by scaling factor for window, see above
    % - divide by number of points in original (NOT zero-padded) data segment
    % - divide by sampling freq to obtain spectral DENSITY
    % - note that cross spectral densities are complex, power spectral densities aren't
    % - the product of a complex number and its conjugate is the same as its absolute
    % value (=magnitude) squared

    % begin by computing power (=auto) spectral densities from each stream
    % type separately:
    % - PSD stream type 1
    disp(['power spectral density (' strmType{1} ')']);
    tmpP1=(tmpDFT1.*conj(tmpDFT1))*scalfac*2/AP.ppSeg*WP.osi*1e-6;
    % average across segments:
    tmpP1=mean(tmpP1,2);
    % - PSD stream type 2
    disp(['power spectral density (' strmType{2} ')']);
    tmpP2=(tmpDFT2.*conj(tmpDFT2))*scalfac*2/AP.ppSeg*WP.osi*1e-6;
    % average across segments:
    tmpP2=mean(tmpP2,2);

    % now iterate thru all channel combinations & compute cross-spec dens
    for ci1=1:AP.nLFPCh
      for ci2=1:AP.nLFPCh
        disp(['cross spectral density (' strmType{1} ' '...
          rawCh(AP.LFPInd(ci1)).nm ', ' strmType{2} ' ' rawCh(AP.LFPInd(ci2)).nm ')']);
        tmpP=(tmpDFT1(:,:,ci1).*conj(tmpDFT2(:,:,ci2)))*scalfac*2/AP.ppSeg*WP.osi*1e-6;
        % average across segments:
        tmpP=mean(tmpP,2);

        % ----- coherence
        tmpCoh=abs(tmpP).^2./(tmpP1(:,ci1).*tmpP2(:,ci2));
        % save unsmoothed
        strm_Coh(:,AP.LFPccInd(ci1),AP.LFPccInd(ci2))=tmpCoh;
        
        % to determine peak, smooth by convoluting with triangular window
        tmpCoh=conv(tmpCoh(pPeakix),krnl);
        % get rid of krnlhPts border points on either end
        tmpCoh=tmpCoh(1+krnlhPts:end-krnlhPts);
        % and since integration requires inthPts on either side of the peaks,
        % look for peaks only in restricted range of tmpCoh
        tmpr=evdeal(tmpCoh(1+inthPts:end-inthPts),'idx','allpeaks');
        % offset-corrected (relative to tmpCoh) positive peak freq
        pp=tmpr.posPeakT{1}+inthPts;
        % if no peak found, entry for that channel will remain nan (default entry of
        % ccDerTemplate)
        if isfinite(pp)
          integr=[];
          for pix=1:length(pp)
            % among the peaks found pick the most prominent one: go a few bins to left
            % and right of the peaks and calculate area
            integr(pix)=sum(tmpCoh([-inthPts:inthPts]+pp(pix)));
          end
          % pick peak with largest integral..
          [m,ix]=max(integr);
          % ..not forgetting to add time offset: pPeakix AND krnlhPts due to convolution
          % (=extension by krnlhPts on either side) & subsequent cutting off of 2*krnlPts
          % on either side
          strm_CohPeakF(1,AP.LFPccInd(ci1),AP.LFPccInd(ci2))=F(pp(ix)+krnlhPts-1+pPeakix(1)-1);
          strm_CohPeak(1,AP.LFPccInd(ci1),AP.LFPccInd(ci2))=tmpr.posPeak{1}(ix);
        end
      end % for ci2=AP.LFPccInd;
    end % for ci1=AP.LFPccInd;
    
    % average coherence in theta freq band
    strm_CohTh=mean(strm_Coh(thFIx,:,:),1);
    
    % now, at the end, give computed results proper name and make fields of r
    tmpSVar=whos('strm_Coh*');
    for tmpVix=1:length(tmpSVar)
      eval(['r(i).' cohType tmpSVar(tmpVix).name(6:end) '=' tmpSVar(tmpVix).name '; clear ' tmpSVar(tmpVix).name]);
    end
  end % if not isempty segmentType
end %for i=1:length(r)

% bring up the summary results figure
drawnow;
