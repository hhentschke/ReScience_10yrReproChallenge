logstr{end+1}='---------- generating audio files..';

% ---- some settings
au_meth='peakdrum';

% resolution of amplitude distribution histograms 
% (needed for proper scaling of the various streams)
res=.0025;
% audio frequency
fs=AP.AudioFs;
% path & file name body of audio file
wfnb=[AP.strmDir '\' DS.abfFn '_'];
      
% ---- preparations
if isstr(AP.rawAudioExcerpt) & strcmpi(AP.rawAudioExcerpt,'full')
  AP.rawAudioExcerpt=[0 discrete2cont(abfi.dataPtsPerChan,osi*1e-6,'intv',1)];
end

% normalize scale
scl=AP.wMagScale/AP.wMagScale(1);
% neuronal signals' sampling frequency
fsNeuro=1e6/osi;
% # of neuronal channels
nAuCh=length(AP.rawChAudioNm);
% indices to neuronal data points 
excPts=cont2discrete(AP.rawAudioExcerpt,osi*1e-6,'intv',1);
% length of audio stream considering fs
% nWavPts=ceil(fs/fsNeuro*nNeuroPts);

% index INTO AP.rawChAnNm to channels to be 'audiolyzed'
for i=1:nAuCh
  auChIdx(i)=strmatch(AP.rawChAudioNm{i},AP.rawChAnNm,'exact');
end


switch au_meth
  % ----- detect theta peaks and replace each peak by a drumbeat; 
  %       gamma and ripples are amplitude-modulated noise/hi freq sounds
  case 'peakdrum'
    % read drum wave & audio sampling freq
    y1=wavread(AP.wffn);
    y1=y1(1:2:end);    
    wc1=round(AP.wCenter/2);
    nBeat1Pts=length(y1);
    % beat #2: same but twice as fast
    y2=y1(1:2:end);
    wc2=round(AP.wCenter/2/2);    
    nBeat2Pts=length(y2);    

    % compute blocks into which data have to be split given the maximum allowed 
    % memory allocation:
    % - each second of neural data requires fsNeuro data points
    % - for each point of neural data, 2*nAuCh*fs/fsNeuro additional points will be generated
    % - each point takes up 8 byte
    % The amount of memory resulting from this should not exceed WP.maxRAM. So,  
    ilen=min(diff(AP.rawAudioExcerpt), 2^20*WP.maxRAM/(8*nAuCh*(fsNeuro+2*fsNeuro*fs/fsNeuro)));
    [intrvls,intrvls_pts]=mkintrvls(AP.rawAudioExcerpt,'resol',osi*1e-6,'ilen',ilen,'olap',0,'border','skip');

    % **** main loop ****    
    for iidx=1:size(intrvls_pts,1)
      % number of pts per block and channel
      nNeuroPts=diff(intrvls_pts(iidx,:))+1;
      % --- gamma:
      % load gamma envelope stream
      for chInd=1:nAuCh
        gaEd(:,chInd)=strmread([AP.strmDir '\' rawCh(auChIdx(chInd)).gammaEnvFn],'intv', intrvls_pts(iidx,:),'verbose',0);  
      end
      % 
      % amplitude distribution for each channel ONLY FROM FIRST DATA BLOCK, assuming that
      % characteristics do not change a lot in the following blocks:
      % disregard amplitudes representing upper .2 % of points since these very
      % likely represent artifacts (which spoil scaling)
      if iidx==1
        gaco=cumh(gaEd,res,'p',[.5 .998]);
      end
      % clip & subtract mean
      for chInd=1:nAuCh
        gaEd(gaEd(:,chInd)>gaco(2,chInd),chInd)=gaco(2,chInd);      
        gaEd(:,chInd)=gaEd(:,chInd)-gaco(1,chInd);
      end
      % interpolate (=upsample to fsAudio)
      for i=nAuCh:-1:1
        w(:,i)=interp1([0:size(gaEd,1)-1]',gaEd(:,i),[(0:nNeuroPts*fs/fsNeuro-1)/(fs/fsNeuro)]');        
      end
      gaEd=[];
      nWavPts=size(w,1);
      % point-by-point multiplication with whatever sound gamma shall represent
      % (channel by channel and even point by point to reduce memory demand)
      % on the occasion, scale 
      for i=nAuCh:-1:1
        for j=1:nWavPts
          w(j,i)=w(j,i)*rand*scl(2);
        end
      end
      
%       % --- ripples:
%       % load ripple envelope stream
%       for chInd=1:nAuCh
%         riEd(:,chInd)=strmread([AP.strmDir '\' rawCh(auChIdx(chInd)).rippleEnvFn],'intv', intrvls_pts(iidx,:),'verbose',0);  
%       end
%       
%       % see comments above
%       if iidx==1
%         % clip all points in between 99.5 % and 99.999 %
%         % (no subtraction of mean here)
%         rico=cumh(riEd,res,'p',[.995 .99999]);
%       end
%       for chInd=1:nAuCh
%         riEd(riEd(:,chInd)<rico(1,chInd),chInd)=0;      
%         riEd(riEd(:,chInd)>rico(2,chInd),chInd)=rico(2,chInd);            
%       end
%       % interpolate (=upsample to fsAudio)
%       for i=nAuCh:-1:1
%         riEdip(:,i)=interp1([0:size(riEd,1)-1]',riEd(:,i),[(0:nNeuroPts*fs/fsNeuro-1)/(fs/fsNeuro)]');
%       end
%       riEd=[];
%       % point-by-point multiplication with whatever sound ripple shall represent
%       for i=nAuCh:-1:1
%         w(1:10:end,i)=w(1:10:end,i)+riEdip(1:10:end,i)*scl(3);
%       end
%       riEdip=[];
      
      % --- theta:
      % load theta stream
      for chInd=1:nAuCh
        thd(:,chInd)=strmread([AP.strmDir '\' rawCh(auChIdx(chInd)).thetaFn],'intv', intrvls_pts(iidx,:),'verbose',0);
      end
      % see comments above
      if iidx==1
        % theta: clip amplitudes representing lower and upper .5 % of values since these very
        % likely represent artifacts (which spoil scaling); subtracting mean not necessary
        thco=cumh(thd,res,'p',[.005 .995]);
      end
      % clip
      for chInd=1:nAuCh
        thd(thd(:,chInd)<thco(1,chInd),chInd)=thco(1,chInd);
        thd(thd(:,chInd)>thco(2,chInd),chInd)=thco(2,chInd);      
      end
      % detect peaks
      tmpr=evdeal(thd,osi,'allpeaks');
      thd=[];
      
      % extend w at both ends to allow for border beats
      w=cat(1,zeros(nBeat1Pts,nAuCh),w,zeros(nBeat1Pts,nAuCh));
      for chInd=1:nAuCh
        % convert output of evdeal to points considering 
        % 1. NEW si 2. the wave's 'acoustic center', and 3. the extended border
        tpos=cont2discrete(tmpr.posPeakT{chInd},1000/fs)-wc2+nBeat1Pts;
        tneg=cont2discrete(tmpr.negPeakT{chInd},1000/fs)-wc1+nBeat1Pts;      
        % for each positive peak implant beat2, for each neagtive one beat1
        for i=1:length(tpos)
          w(tpos(i):tpos(i)+nBeat2Pts-1,chInd)=w(tpos(i):tpos(i)+nBeat2Pts-1,chInd)+y2*tmpr.posPeak{chInd}(i);
        end
        for i=1:length(tneg)
          w(tneg(i):tneg(i)+nBeat1Pts-1,chInd)=w(tneg(i):tneg(i)+nBeat1Pts-1,chInd)+y1* -tmpr.negPeak{chInd}(i);
        end
      end
      % cut down again
      w([1:nBeat1Pts end-nBeat1Pts+1:end],:)=[];
      
      % --- concluding section:  
      % normalize wav (assuming there is no offset)
      ma=max(w,[],1);    mi=min(w,[],1);
      % normalize such that max peak corresponds to +/-.99
      mm=max(abs([ma; mi]),[],1);
      w=.99*w./repmat(mm,size(w,1),1);
      wfn=[wfnb sprintf('%.2i',iidx) '.wav'];
      wavwrite(w,fs,16,wfn);
      wavplay(w,fs,'async');
      w=[];
    end
    
    
  case 'awfm' 
    % ----- arbitrary waveform frequency modulation
    % read basic wave & audio sampling freq
    [y,fs]=wavread(AP.wffn);
    nWavPts=length(y);
    
    
    % regard the lower and upper .375 % as outliers
    [nix,i]=min(abs(cs-.005));
    rng(1)=ab(min(i));
    [nix,i]=min(abs(cs-.995));
    rng(2)=ab(max(i));
    % scale: 99.x % of signal should be bound by [1-AP.modDepth 1]
    thd=(thd-rng(1))/diff(rng)*AP.modDepth+(1-AP.modDepth);
    
    
    % interpolation of neural data (=upsampling to fsAudio
    for i=1:size(tmpd,2)
      tmpd_ip(:,i)=interp1([0:size(tmpd,1)-1]',tmpd(:,i),[(0:(nNeuroPts-1)*fs/fsNeuro)/(fs/fsNeuro)]');
    end
    clear tmpd;

%     % interpolation of neural data (=upsampling to fsAudio)
%     for i=1:nAuCh
%       thd_ip(:,i)=interp1([0:nNeuroPts-1]',thd(:,i),[(0:(nNeuroPts-1)*fs/fsNeuro)/(fs/fsNeuro)]');
%     end
%     clear thd;
%     nWavPts=size(thd_ip,1);
    
    
    w=awfm(tmpd_ip,fs,y);
    
    clear tmp* ah ab cs nix
end