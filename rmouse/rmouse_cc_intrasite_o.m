% code below deletes faulty elements resulting from non-preallocation of
% results fields thgaeCCZScore and thgaeCCZTestP  up until dec 2005

% for i=1:length(r)
%   if ~isempty(r(i).iPts)
%     disp([r(i).segmentType '..']);
%     if size(r(i).thgaeCCZScore,1)>r(i).ni
%       disp('deleting faulty entries')
%       r(i).thgaeCCZScore(r(i).ni:end,:)=[];
%       r(i).thgaeCCZTestP(r(i).ni:end,:)=[];
%     end
%   end
% end
      

% catch some silly errors
if isempty(AP.ccNShuffle)
  warning('AP.ccNShuffle is empty - set to 0 to prevent shuffling');
end

pack
switch ccType
  case 'thgae'
    strmType={'theta','gammaEnv'};
    sChInd=AP.OPT_thgaeCCStartChInd;
    attractLag=AP.OPT_thgaeCCAttractLag;
    % detect peaks within this interval centered around t=0
    ccw=cont2discrete(AP.thccw,osi*.001,'intv',0)-1; % sic!
    cca=AP.cca;
    % at the fissure, where theta and gammaEnv are most strongly correlated, the
    % correlation peak is negative. Therefore, the algorithm should hop into the
    % trough
    ccPolarity='negative';
  case 'dethe'
    strmType={'delta','thetaEnv'};    
    sChInd=AP.OPT_detheCCStartChInd;
    attractLag=AP.OPT_detheCCAttractLag;
    % detect peaks within this interval centered around t=0
    ccw=cont2discrete(AP.deccw,osi*.001,'intv',0)-1; % sic!
    cca=AP.deccw;    
    % deal with positive cc peaks
    ccPolarity='positive';
  otherwise
    error('illegal ccType in intra-site CC computation');
end

% catch typos and set index for temporary vars to positive and negative peaks
switch ccPolarity
  case 'negative'
    tmpPeakIx=2;
  case 'positive'
    tmpPeakIx=1;
  otherwise
    error('illegal choice for ccPolarity');
end


if isempty(sChInd)  
  sChInd=AP.LFPpcInd1;  
end

% generate order of channels to analyze: start with channel determined above, 
% work way down, then up
if nLFPCh>1
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
% tell me what fraction of the CC output (length 2*AP.ccLagPts+1) the excerpts
% for peak & lat calc take up (needed for guesses of memory demand)
tmpSFrac=diff(ccw)/(2*AP.ccLagPts+1);
for i=1:length(r)
  if ~isempty(r(i).iPts)
    disp([r(i).segmentType '..']);
    % preallocate results variables  
    % I. 2D array, each column the mean CC derived from one pair of streams of one channel 
    cccMn=repmat(nan,[2*AP.ccLagPts+1 nLFPCh]);
    cccStd=cccMn;
    % II. amplitudes and lags of designated (single) peak in individual CC segments.
    cccPeak=repmat(nan,[r(i).ni  nLFPCh]);
    cccPeakT=repmat(nan,[r(i).ni  nLFPCh]);
    % III. peak decay and Z score
    if strcmpi(ccType,'thgae')
      cccPosPeakDecay=repmat(nan,[r(i).ni  nLFPCh]);
      cccPosPeakDecayMn=repmat(nan,[1  nLFPCh]);
      cccPosPeakDecayStd=repmat(nan,[1  nLFPCh]);
      cccNegPeakDecay=cccPosPeakDecay;
      cccNegPeakDecayMn=repmat(nan,[1  nLFPCh]);
      cccNegPeakDecayStd=repmat(nan,[1  nLFPCh]);
      cccZScore=cccPosPeakDecay;
      cccZTestP=cccPosPeakDecay;
    end
    % intermediate results var: for each channel, CC of individual segments in columns
    segCC=repmat(nan,[2*AP.ccLagPts+1 r(i).ni]);
    % intermediate: cell arrays holding peaks of current channel;
    % 2D (no shuffling) or 3D; 1st row positive-going peaks, second row
    % negative-going
    tmpPeak=cell([2 r(i).ni AP.ccNShuffle+1]);
    tmpPeakT=tmpPeak;          
    % assess size taken up by major variables considering/assuming the following:
    % - one empty cell of a cell array takes up 40=5*8 bytes
    % - number of peaks in the CC = number expected from mean theta freq
    % - number of peaks in the CC from shuffled data = twice the above
    % size of one cell containing peak amplitudes of one individual data
    % segment, expressed in double scalar equivalents
    tmpsz=(5+mean(AP.theta)*discrete2cont(length(cci),osi*1e-6));
    tmpVarSz=...
      r(i).ni*2*2*2*tmpsz+...        % tmpr (no shuffling, all segs at once) & tmpPeak, tmpPeakT which together are the size of tmpr
      r(i).ni*2*nLFPCh+...             % .thgaeCCPeak & .thgaeCCPeakT
      2*prod(size(cccMn))+... % mean & std
      prod(size(segCC))*(1+2*tmpSFrac);  % last factor accounts for mem demand of evdeal
    if AP.ccNShuffle          
      % intermediate results var: peak amplitudes of CC from all shuffled data
      tmpShPeak=repmat(nan,[r(i).ni  1  AP.ccNShuffle]);
      % intermediate results var: CC between theta and a) ONE segment of gamma 
      % envelope (in 1st column) and b) phase-shuffled versions of it (2nd and up
      % columns)
      sh_segCC=repmat(nan,[2*AP.ccLagPts+1 1 AP.ccNShuffle+1]);
      % assess size taken up by variables:
      tmpVarSz=...
        2*2*(2*AP.ccNShuffle+1)*tmpsz+...    % tmpr (shuffling, 2*AP.ccNShuffle because more peaks expected in shuffled data) 
        r(i).ni*2*tmpsz+...                  % tmpPeak, tmpPeakT
        r(i).ni*2*nLFPCh+...                 % .thgaeCCPeak & .thgaeCCPeakT
        r(i).ni*(2*AP.ccNShuffle+1)+...      % tmpShPeak
        2*prod(size(cccMn))+...     % mean & std
        prod(size(segCC))*(1+2*tmpSFrac)+...   % last factor accounts for mem demand of evdeal                         
        prod(size(sh_segCC))*(1+2*tmpSFrac);    % last factor accounts for mem demand of evdeal
      % prepare temporary intermediate results fig
      % **** plot suspended April 06
      if 0
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
    end            
    disp([ 'peak additional memory demand ~ ' num2str(tmpVarSz*8/2^20,'%4.1f') ' Mb']);
    % reset list index and attractT (expected location of peaks) 
    li=0;
    if isempty(attractLag)    
      attractT=nan*(1:nLFPCh);
    end
    for chInd=chList
      % chInd is the index into all results variables whereas AP.LFPInd(chInd) is the
      % index to use with rawCh
      li=li+1;
      disp([r(i).segmentType ': XC ' rawCh(AP.LFPInd(chInd)).nm ' ' strmType{1} ' vs ' strmType{2} ' ..']);
      % load data
      eval(['tmpd1=strmread([AP.strmDir ''\'' rawCh(AP.LFPInd(chInd)).' strmType{1} 'Fn],''nPts'',r(i).iPts(end,2),''verbose'',0);']);
      eval(['tmpd2=strmread([AP.strmDir ''\'' rawCh(AP.LFPInd(chInd)).' strmType{2} 'Fn],''nPts'',r(i).iPts(end,2),''verbose'',0);']);      
      if AP.ccNShuffle
        disp('computing CC from original and shuffled segments..')
      else
        disp('computing CC from segments..')
      end
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
        end   % if:shuffle  
      end % for:segments

%       % uncomment block below below for export of raw CC
%       tmp_time=floor(mean(r(i).iPts,2));
%       save([AP.resPath '\' ccType '_' AP.resFn '_' rawCh(AP.LFPInd(chInd)).nm '_' r(i).segmentType],'segCC','tmp_time','osi','abfi','AP');

      disp('picking peaks..')
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

      % alas, we need a loop to determine peak decay
      if strcmpi(ccType,'thgae')
        for cpi=1:r(i).ni
          % determine peak decay: sort according to amplitude and divide second
          % largest by largest
          tmpSortP=sort(tmpPeak{1,cpi,1});
          if length(tmpSortP)>=2 & all(isfinite(tmpSortP))
            cccPosPeakDecay(cpi,chInd)=tmpSortP(end-1)/tmpSortP(end);
          end
          tmpSortP=sort(tmpPeak{2,cpi,1});
          if length(tmpSortP)>=2 & all(isfinite(tmpSortP))
            cccNegPeakDecay(cpi,chInd)=tmpSortP(2)/tmpSortP(1);
          end
        end
        OKix=isfinite(cccPosPeakDecay(:,chInd));
        cccPosPeakDecayMn(1,chInd)=mean(cccPosPeakDecay(OKix,chInd));
        cccPosPeakDecayStd(1,chInd)=std(cccPosPeakDecay(OKix,chInd));
        OKix=isfinite(cccNegPeakDecay(:,chInd));
        cccNegPeakDecayMn(1,chInd)=mean(cccNegPeakDecay(OKix,chInd));
        cccNegPeakDecayStd(1,chInd)=std(cccNegPeakDecay(OKix,chInd));
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
          % since there may be original or shuffled segments without a peak
          % in the specified interval, there may be nans, so this is a good
          % place to compute the z score/test
          ix1=isfinite(tmpShPeak(k,1,:));
          if isfinite(cccPeak(k,chInd)) & sum(ix1)>1
            [nada,p,nada2,zval]=ztest(cccPeak(k,chInd),mean(tmpShPeak(k,1,ix1)),std(tmpShPeak(k,1,ix1)));
          else
            zval=nan;
            p=nan;
          end
          cccZScore(k,chInd)=zval;
          cccZTestP(k,chInd)=p;
          % cccZScore(k,chInd)=(cccPeak(k,chInd)-mean(tmpShPeak(k,1,ix1)))/std(tmpShPeak(k,1,ix1));
        end % if:shuffle
      end  % for:segments    
      % **** plot suspended April 06
      if 0 && AP.ccNShuffle 
        % plot 
        [co,h,bins]=cumh(shiftdim(tmpShPeak(:,1,:),2),.01);
        subplot(8,ceil(nLFPCh/8),chInd);
        plot(bins,h);
        hold on;
        % replicating & arranging data as columns results in individual points 
        % being plotted in colors matching those of cum histograms
        plot(repmat((cccPeak(:,chInd))',2,1),ones(2,r(i).ni)*1.08,'v');
        nicexax;
        if chInd==1, title([r(i).segmentType ', ' rawCh(AP.LFPInd(chInd)).nm]); 
        else title(rawCh(AP.LFPInd(chInd)).nm); 
        end
        drawnow
      end % if:shuffle
      % do all these calculations channel by channel - easier to kick out nans
      cccPeakMn(1,chInd)=mean(cccPeak(OKix,chInd),1);
      cccPeakStd(1,chInd)=std(cccPeak(OKix,chInd),0,1);                
      cccPeakTMn(1,chInd)=mean(cccPeakT(OKix,chInd),1);
      cccPeakTStd(1,chInd)=std(cccPeakT(OKix,chInd),0,1);
    end % for:channels
    eval(['r(i).' ccType 'CCMn=cccMn;']);
    eval(['r(i).' ccType 'CCStd=cccStd;']);
    eval(['r(i).' ccType 'CCPeak=cccPeak;']);
    eval(['r(i).' ccType 'CCPeakT=cccPeakT;']);
    if exist('cccZScore','var')
      eval(['r(i).' ccType 'CCZScore=cccZScore;']);
      eval(['r(i).' ccType 'CCZTestP=cccZTestP;']);
    end
    eval(['r(i).' ccType 'CCPeakMn=cccPeakMn;']);
    eval(['r(i).' ccType 'CCPeakStd=cccPeakStd;']);    
    eval(['r(i).' ccType 'CCPeakTMn=cccPeakTMn;']);
    eval(['r(i).' ccType 'CCPeakTStd=cccPeakTStd;']);  
    if strcmpi(ccType,'thgae')
      eval(['r(i).' ccType 'CCPosPeakDecay=cccPosPeakDecay;']);
      eval(['r(i).' ccType 'CCPosPeakDecayMn=cccPosPeakDecayMn;']);
      eval(['r(i).' ccType 'CCPosPeakDecayStd=cccPosPeakDecayStd;']);
      eval(['r(i).' ccType 'CCNegPeakDecay=cccNegPeakDecay;']);
      eval(['r(i).' ccType 'CCNegPeakDecayMn=cccNegPeakDecayMn;']);
      eval(['r(i).' ccType 'CCNegPeakDecayStd=cccNegPeakDecayStd;']);
    end
  end % if: ~isempty(segType)
end % for:segmentType

% **** plot suspended April 06
% if AP.ccNShuffle, close(fhShCC); end

clear ccw cci ccc* segCC sh_segCC tmp* fhShCC i ii k ix* n m h co bins li cmpList chList p1 p2 p rV sChInd OKix
pack