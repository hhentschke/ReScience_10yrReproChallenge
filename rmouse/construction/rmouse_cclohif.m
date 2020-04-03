% detect peaks within this interval centered around t=0
ccw=cont2discrete(AP.thccw,osi*.001,'intv',0)-1; % sic!
cca=AP.cca;
cci=AP.ccLagPts+[ccw(1):ccw(2)];

% the max amount of data (in seconds) that should be used for this (extremely
% time-consuming) analysis
maxAnT=200;

% at the fissure, where theta and gammaEnv are most strongly correlated, the
% correlation peak is negative. Therefore, the algorithm should hop into the
% trough
ccPolarity='negative';
tmpPeakIx=2;

% lo freq in steps of 1 Hz, span 4 Hz
fr1=[2:11]';
fr1(:,2)=fr1(:,1)+4;

% hi freq in steps of 5 Hz, span 15 Hz
fr2=[5:5:160]';
fr2(:,2)=fr2(:,1)+15;

% half the original sampling rate in Hz
hfsample=.5*1e6/si;
% filter order (=steepness)
ord=3;

% order of channels: first s.p. (located ~500 um above s.lm., then s.lm.)
chList=AP.LFPpcInd1+[-5 0];

fh_ccLoHiF=mkfig('ccLoHiF');
clf, orient portrait
labelscale('fontSz',8,'scaleFac',1.0,'lineW',1.0,'markSz',4);

for i=1:length(r)
  % exploring only
  if ~isempty(r(i).iPts) & i==explV
    disp([r(i).segmentType '..']);
    % number of segments to be analyzed here
    nAnSeg=maxAnT/(AP.ppSeg*osi*1e-6);
    % indices into these
    if nAnSeg<r(i).ni
      segIx=ceil(r(i).ni:-r(i).ni/nAnSeg:1);
    else
      segIx=r(i).ni:-1:1;
    end
    
    % preallocate results variables  
    % I. 2D array, each column the mean CC derived from one pair of streams of one channel 
    cccMn=repmat(nan,[2*AP.ccLagPts+1 length(chList)]);
    cccStd=cccMn;
    for chInd=chList
      % chInd is the index into all results variables whereas AP.LFPInd(chInd) is the
      % index to use with rawCh
      disp([r(i).segmentType ': XC ' rawCh(AP.LFPInd(chInd)).nm ' LO vs HI Env ..']);
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
      for frix1=length(fr1):-1:1
        disp(['lo: ' num2str(fr1(frix1,:)) ' Hz']);
        % edge frequencies of passband normalized to hfsample, the Nyquist freq
        wn=fr1(frix1,:)/hfsample;
        [b_lo,a_lo]=butter(ord,wn);
        for frix2=length(fr2):-1:1
          disp(['hi: ' num2str(fr2(frix2,:)) ' Hz']); 
          % edge frequencies of passband normalized to hfsample, the Nyquist freq
          wn=fr2(frix2,:)/hfsample;
          [b_hi,a_hi]=butter(ord,wn);
          % loop only over a subset of segments (the equivalent of 3 min)
          for k=segIx
            % indices into both streams for current segment
            tmpIdx=r(i).iPts(k,1):r(i).iPts(k,2);
            % filter (& envelope) each segment individually (which is faster if the
            % total time covered by the segments is much smaller than the whole
            % stream)
            loD=filtfilt(b_lo,a_lo,tmpd(tmpIdx));
            hiD=filtfilt(b_hi,a_hi,tmpd(tmpIdx));
            hiD=detrend(abs(hilbert(hiD)),'constant');
            % CC
            segCC(:,k)=xxcorr(loD,hiD,AP.ccLagPts,AP.ccScaleOpt);
          end % for:segments
          % ------- mean & std of CC waveform   **segCC is not emptied**
          cccMn(:,chInd)=mean(segCC,2);
          % cccStd(:,chInd)=std(segCC,0,2);
          tmpr=evdeal(cccMn(cci,chInd),osi,'minmaxpeak');
          
          if strcmpi(ccPolarity,'negative')
            ccm(frix1,frix2)= -1*tmpr.minPeak;
            % peakT=tmpr.negPeakT{1}+discrete2cont(ccw(1),osi*.001,'intv',1);
          else 
            error('sdf');
          end
        end % for:hi freqs
      end% for:lo freqs
      % plot!
      subplot(length(chList),1,find(chList==chInd)),
      imagesc(mean(fr2,2),mean(fr1,2),ccm,[.0 .05]);
      xlabel('HI Frequency (Hz)');
      ylabel('LO Frequency (Hz)');      
      colorbar;
      drawnow;
      saveas(gcf,[figName '_ccLoHiF'],'fig');
    end % for:chInd
  end % if: ~isempty(segType)
end % for:segmentType


clear ccw cci ccc* segCC shnix tmp* fhShCC i ii k ix* n m h co bins li cmpList chList p1 p2 p rV sChInd OKix
pack