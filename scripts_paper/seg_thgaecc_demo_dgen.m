% This routine generates data used for plots illustrating details of
% crosscorrelation theta-gammaEnv & Matt's shuffle statistics. It is
% directly derived from rmouse_ccIntraSite. *** must be called from within
% rmouse as a special job *** *** make sure that in the calling AP only one
% channel is requested for analysis ***

% --- specific settings (i) apart from those specified in, or (ii) overriding those in AP
behav={'exploring'};
% excerpt of data to display in raw fig
intv=[99.55 100.8]+45; % a realistic picture: similarity, but not perfect

% --- general preps
strmType={'theta','gammaEnv'};
sChInd=AP.LFPpcInd1;

attractLag=AP.OPT_thgaeCCAttractLag;
% detect peaks within this interval centered around t=0
ccw=cont2discrete(AP.thccw,osi*.001,'intv',0)-1; % sic!
cca=AP.cca;
% at the fissure, where theta and gammaEnv are most strongly correlated, the
% correlation peak is negative. Therefore, the algorithm should hop into the
% trough
ccPolarity='negative';
tmpPeakIx=2;

% generate order of channels to analyze: start with channel determined above, 
% work way down, then up
if nLFPCh>1
  error([mfilename ' makes sense only with one channel to be analyzed']);
  % order of channels 
  chList=[sChInd:-1:1, sChInd+1:nLFPCh];
  % order of channels to compare with:
  cmpList=[nan, sChInd:-1:2, sChInd:nLFPCh-1];
else 
  chList=sChInd;
  cmpList=nan;
end

% don't change - original values of attractLag (that is, the question whether it was empty
% or not) is needed below
if ~isempty(attractLag)
  attractT=attractLag;
end

cci=AP.ccLagPts+[ccw(1):ccw(2)];

% ** exploring only **
for i=explV
  if ~isempty(r(i).iPts)
    disp([r(i).segmentType '..']);
    % preallocate results variables  
    % I. 2D array, each column the mean CC derived from one pair of streams of one channel 
    cccMn=repmat(nan,[2*AP.ccLagPts+1 nLFPCh]);
    cccStd=cccMn;
    % II. amplitudes and lags of designated (single) peak in individual CC segments.
    cccPeak=repmat(nan,[r(i).ni  nLFPCh]);
    cccPeakT=repmat(nan,[r(i).ni  nLFPCh]);
    % intermediate results var: for each channel, CC of individual segments in columns
    segCC=repmat(nan,[2*AP.ccLagPts+1 r(i).ni]);
    % intermediate: cell arrays holding peaks of current channel;
    % 2D (no shuffling) or 3D; 1st row positive-going peaks, second row
    % negative-going
    tmpPeak=cell([2 r(i).ni AP.ccNShuffle+1]);
    tmpPeakT=tmpPeak;          
    if AP.ccNShuffle          
      % intermediate results var: peak amplitudes of CC from all shuffled data
      tmpShPeak=repmat(nan,[r(i).ni  1  AP.ccNShuffle]);
      % intermediate results var: CC between theta and a) ONE segment of gamma 
      % envelope (in 1st column) and b) phase-shuffled versions of it (2nd and up
      % columns)
      sh_segCC=repmat(nan,[2*AP.ccLagPts+1 1 AP.ccNShuffle+1]);
      % prepare temporary intermediate results fig
      ftag='XC intra-site: significance of peaks';
      fhShCC=findobj('tag',ftag);
      if isempty(fhShCC), fhShCC=figure;
      else  figure(fhShCC);
      end
      tmpScrSz=get(0,'Screensize');
      tmpScrSz(1)=tmpScrSz(1)+diff(tmpScrSz([1 3]))*.5;
      tmpScrSz(2)=tmpScrSz(2)+diff(tmpScrSz([2 4]))*.05;  
      tmpScrSz(3)=diff(tmpScrSz([1 3]));              
      tmpScrSz(4)=diff(tmpScrSz([2 4]))*.9;                
      set(fhShCC,'position',tmpScrSz,'tag',ftag,'name',ftag,...
        'color',[0.87 0.87 0.87],'numbertitle','off');
      clf;
    end            
    % reset list index and attractT (expected location of peaks) 
    li=0;
    if isempty(attractLag)    
      attractT=nan*(1:nLFPCh);
    end
    % identify segment containing data excerpt to be plotted. this is also
    % the segment whose shuffled versions will be displayed 
    tmpt=intv(1)*1000-discrete2cont(r(i).iPts(:,1),osi*.001,'intv',0);
    [nada,pickSegIx]=min((tmpt>=0).*tmpt);
    pickSegIx=pickSegIx-1;
    
    for chInd=chList
      % chInd is the index into all results variables whereas AP.LFPInd(chInd) is the
      % index to use with rawCh
      li=li+1;
      disp([r(i).segmentType ': XC ' rawCh(AP.LFPInd(chInd)).nm ' ' strmType{1} ' vs ' strmType{2} ' ..']);
      % load data
      eval(['tmpd1=strmread([AP.strmDir ''\'' rawCh(AP.LFPInd(chInd)).' strmType{1} 'Fn],''nPts'',r(i).iPts(end,2),''verbose'',0);']);
      eval(['tmpd2=strmread([AP.strmDir ''\'' rawCh(AP.LFPInd(chInd)).' strmType{2} 'Fn],''nPts'',r(i).iPts(end,2),''verbose'',0);']);
      for k=r(i).ni:-1:1
        % indices into both streams for current segment
        tmpIdx=r(i).iPts(k,1):r(i).iPts(k,2);
        % detrended theta
        tmpd1_exc=detrend(tmpd1(tmpIdx),'constant');
        % gamma envelope
        tmpd2_exc=tmpd2(tmpIdx);
        % CC 'real' gamma envelope
        segCC(:,k)=xxcorr(tmpd1_exc,detrend(tmpd2_exc,'constant'),AP.ccLagPts,AP.ccScaleOpt);
        if AP.ccNShuffle
          % sorted amplitude of original
          tmp_sortAmp=sort(tmpd2_exc);
          % DFT of original
          tmp_offt=fft(tmpd2_exc);
          % generate phase shuffled versions of envelope:
          tmp4=shuffle(1,repmat(tmpd2_exc,1,AP.ccNShuffle));
          % convert from magnitude/phase to real/imaginary part
          [tmp_re,tmp_im]=pol2cart(angle(fft(tmp4)),repmat(abs(tmp_offt),1,AP.ccNShuffle));
          tmp3=tmp_re+sqrt(-1)*tmp_im;
          % inverse DFT: tmp_shEnv
          tmp_shEnv=real(ifft(tmp3));
          % final step in surrogate data generation + CC: need a loop now
          for ii=1:AP.ccNShuffle
            % sort according to amplitude & replace amplitudes by original
            % segment's amplitude
            tmp8=sortrows([tmp_shEnv(:,ii) [1:AP.ppSeg]'],1);
            tmp_shEnv(tmp8(:,2),ii)=tmp_sortAmp;
            sh_segCC(:,1,ii+1)=xxcorr(tmpd1_exc,detrend(tmp_shEnv(:,ii),'constant'),AP.ccLagPts,AP.ccScaleOpt);                
          end
          % put 'real' CC as first column
          sh_segCC(:,1,1)=segCC(:,k);
          % at this point, detect and keep peaks of CC from real and surrogate data - 
          % the relevant peak(s) will be fished out below
          tmpr=evdeal(sh_segCC(cci,1,:),osi,'allpeaks');
          % ---------- pick positive-going peaks from real and shuffled pairs
          tmpPeak(1,k,:)=tmpr.posPeak;
          tmpPeakT(1,k,:)=tmpr.posPeakT;
          % ---------- same thing for negative-going peaks 
          tmpPeak(2,k,:)=tmpr.negPeak;
          tmpPeakT(2,k,:)=tmpr.negPeakT;
          
          if k==pickSegIx
            % save 
            thetaD=tmpd1_exc;
            gammaEnvD=tmpd2_exc;
            shGammaEnvD=tmp_shEnv;
            si=osi;
            save([AP.resPath '\' mfilename '.mat'],'si','thetaD','gammaEnvD','shGammaEnvD','sh_segCC');
          end
          
        end   % if:shuffle  
      end % for:segments
      if ~AP.ccNShuffle
        tmpr=evdeal(segCC(cci,:),osi,'allpeaks');
        % ---------- pick positive-going peaks 
        tmpPeak(1,:)=tmpr.posPeak;
        tmpPeakT(1,:)=tmpr.posPeakT;
        % ---------- same thing for negative-going peaks 
        tmpPeak(2,:)=tmpr.negPeak;
        tmpPeakT(2,:)=tmpr.negPeakT;
      end % if:~shuffle
      % ------- mean & std of CC waveform   **segCC is not emptied**
      cccMn(:,chInd)=mean(segCC,2);
      cccStd(:,chInd)=std(segCC,0,2);      
      % based on comparisons of crosscorrelations of neighboring pairs of
      % theta & gammaEnv streams, determine where within original cc interval
      % to fish for peaks
      tmpr=evdeal(cccMn(cci,chInd),osi,'allpeaks');
      
      if strcmpi(ccPolarity,'negative')
        peak= -1*tmpr.negPeak{1};
        peakT=tmpr.negPeakT{1}+discrete2cont(ccw(1),osi*.001,'intv',1);
      else
        peak=tmpr.posPeak{1};
        peakT=tmpr.posPeakT{1}+discrete2cont(ccw(1),osi*.001,'intv',1);
      end

      if li==1
        % first run: the relevant peak (=the largest) should be unambiguous. 
        % So, in mean CC waveform find T of that peak within 
        % designated interval..
        ix1=find(peakT>=cca(1) & peakT<=cca(2));
        [m,ix2]=max(peak(ix1));
        % the location of this peak is the attractor for peaks to fish out in
        % individual segments
        if isempty(attractLag)
          attractT(chInd)=peakT(ix1(ix2));        
        end
      else
        % in pairs of cc functions, find peak in 2nd cc closest to 'attractor' peak 
        % in 1st (not heeding sign of lag or amplitude). The fact that the closest 
        % peak is looked for implies that we dont accept phase shifts of more than 
        % half a typical period
        if isempty(attractLag)
          [m,ix1]=min(abs(attractT(cmpList(li))-peakT));
          attractT(chInd)=peakT(ix1);
        end
      end
      % local, channel specific cca (shifts with attracting peak)
      shcca=cca+attractT(chInd);
      % indices into temp variables with good (=non-nan) values for lag and
      % amplitude of peaks
      OKix=[];
      for k=r(i).ni:-1:1
        for shi=1:size(tmpPeak,3)
          % 1. local vars & subtract offset  
          % ** note that here the choice between pos and neg peaks lies in var
          % tmpPeakIx
          peak=abs(tmpPeak{tmpPeakIx,k,shi});
          peakT=tmpPeakT{tmpPeakIx,k,shi}+discrete2cont(ccw(1),osi*.001,'intv',1);              
          if isempty(peak)
            if shi>1
              estr=['channel ' rawCh(AP.LFPInd(chInd)).nm ', segment # ' int2str(k),...
                  ', shuffled version # ' int2str(shi),... 
                  ': CC ' strmType{1} '-' strmType{2} ' does not have a single peak within [',...
                  int2str(AP.ccLagPts*[-1 1]) '] ticks'];
            else
              estr=['channel ' rawCh(AP.LFPInd(chInd)).nm ', segment # ' int2str(k),...
                  ': CC ' strmType{1} '-' strmType{2} ' does not have a single peak within [',...
                  int2str(AP.ccLagPts*[-1 1]) '] ticks'];
            end
            error(estr);
          end
          ix1=find(peakT>=shcca(1) & peakT<=shcca(2));       
          % unlikely, but possible: no peak found
          if isempty(ix1)
            % if we're dealing with shuffled data..
            % ** note the -1 in index: first layer of tmpPeak = unshuffled values 
            if shi>1
              tmpShPeak(k,1,shi-1)=nan;
            else
              cccPeakT(k,chInd)=nan;
              cccPeak(k,chInd)=nan;
            end
          else %if:<no peak found>
            [m,ix2]=min(abs(peakT(ix1)-attractT(chInd)));                  
            if shi>1
              % just retain amplitude
              tmpShPeak(k,1,shi-1)=peak(ix1(ix2));
            else
              OKix=[OKix k];
              cccPeakT(k,chInd)=peakT(ix1(ix2));
              cccPeak(k,chInd)=peak(ix1(ix2));
            end
          end % if:<no peak found>  
        end  % for:shuffled segments
        if AP.ccNShuffle 
          % since there may be shuffled segments without a peak in the specified
          % interval, there may be nans, so this is a good place to compute the 
          % z score/test
          ix1=isfinite(tmpShPeak(k,1,:));
          [foo,p,ci,zval]=ztest(cccPeak(k,chInd),mean(tmpShPeak(k,1,ix1)),std(tmpShPeak(k,1,ix1)));
          cccZScore(k,chInd)=zval;
          cccZTestP(k,chInd)=p;
        end % if:shuffle
      end  % for:segments    
      
      % cum histograms of CC peak values for each of the shuffled segments
      [co,shCumH,shBins]=cumh(shiftdim(tmpShPeak(:,1,:),2),.005);
      [co,realCumH,realBins]=cumh(cccPeak(:,chInd),.005);
      % save
      save([AP.resPath '\' mfilename '.mat'],'cccZScore','cccZTestP','shCumH','shBins','realCumH','realBins','-append');
    end % for:channels
  end % if: ~isempty(segType)
end % for:segmentType
if AP.ccNShuffle, close(fhShCC); end


clear ccw cci ccc* segCC sh_segCC tmp* fhShCC i ii k ix* n m h co bins li cmpList chList p1 p2 p rV sChInd OKix
pack