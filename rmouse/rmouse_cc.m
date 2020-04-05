function rmouse_cc(rawCh,strmType,STshort)
% Strategy for determination of peaks and their lags:
% - compute average of segment-wise CC between first channel and its nearest 
%   neighbors
% - identify POI (peak of interest) in this CC: the one with the largest amplitude
% - the lag of this peak makes an entry in matrix expLag (holding expected lags)
% - then, in each individual CC (resulting from the segment-wise computations) find
%   the peak closest to POI; collect all peaks and compute their mean lag. This mean 
%   lag replaces the previous entry in expLag
% - repeat this for all pairs of immediately adjacent channels
% - deal with pairs of not immediately adjacent channels: the lag between them should
%   be the sum of all imm adj channel pairs connecting them. E.g., the lag between
%   ch3 and ch5 should be the sum of lags ch3-ch4 and ch4 ch5. So, compute the expected 
%   lag, and pick the peaks in the segment-wise CC which come closest to it (regardless 
%   of amplitude)
% This strategy generally works well; it gets derailed if there are large phase jumps
% between adjacent channels or if the true peaks are buried within multiple other peaks 
% due to line hum: in this case the algorithm will eventually jump on a line hum peak
% (the closest it'll get), offsetting the expected lag and transmitting this error.

global DS AP WP r

warning('off','stats:regress:RankDefDesignMat');

switch strmType
  case 'delta'
    % detect peaks within this interval centered around t=0
    ccw=cont2discrete(AP.deccw,WP.osi*.001,'intv',0)-1; % sic!
  case 'theta'
    ccw=cont2discrete(AP.thccw,WP.osi*.001,'intv',0)-1; % sic!
  case 'thetaHiEnv'
    ccw=cont2discrete(AP.deccw,WP.osi*.001,'intv',0)-1; % sic!
  case 'thetaLoEnv'
    ccw=cont2discrete(AP.ccLagPts*[-.95 .95],WP.osi*.001,'intv',0)-1; % sic!
  case {'gamma','gammaEnv'}
    ccw=cont2discrete(AP.gaccw,WP.osi*.001,'intv',0)-1; % sic!
  case {'gammaNarrow','gammaNarrowEnv'}
    ccw=cont2discrete(AP.gaNaccw,WP.osi*.001,'intv',0)-1; % sic!
otherwise
    error('illegal streamType in CC computation');
end

% --- some parameters
maxFracd=.4;

% --- interim figure 
labelscale('fontSz',10,'scaleFac',1.0,'lineW',1.0,'markSz',5); 
tmpftag='CCLags';
fh_ccLag=findobj('tag',tmpftag);
if isempty(fh_ccLag), fh_ccLag=figure;
else  figure(fh_ccLag);
end
tmpScrSz=get(0,'Screensize');
tmpScrSz=round([tmpScrSz(3)*.05  tmpScrSz(4)*.45  tmpScrSz(3)*.9  tmpScrSz(4)*.4]);  
set(fh_ccLag,'position',tmpScrSz,'tag',tmpftag,'name','peak detection in CC' ,...
  'color',[0.8 0.8 0.8],'numbertitle','off','menubar','none');

% --- prelims
cci=AP.ccLagPts+[ccw(1):ccw(2)];
% accounts for lag offset of cc (in computation of peak occurrence times in cc functions)
tOffs=discrete2cont(ccw(1),WP.osi*.001,'intv',0);
% the number of diagonals in cc matrix
nDiags=length(WP.diagNccix);

% --- load it up, Scotty (if you can)
% if RAM permits, upload all channels completely here (thus circumventing multiple 
% read operations with each segment type and channel)
mame=max(cat(1,r(:).dmem));
mani=max(cat(1,r(:).lastPt));
tmpD=[];  

if str2num(WP.mver)>=7.6
  userview=memory;
  isMemAvailable=userview.MaxPossibleArrayBytes/2^20 > mame*(AP.nLFPCh+2)/AP.nLFPCh;
else
  if str2num(WP.mver)>=7 && WP.javaEnabled
    feature('MemStats');
    isMemAvailable=feature('DumpMem')/2^20 > mame*(AP.nLFPCh+2)/AP.nLFPCh;
  else
    isMemAvailable=logical(1);
    try
      tmpD=repmat(0,[mani AP.nLFPCh]);
    catch
      isMemAvailable=logical(0);
    end
  end
end

if isMemAvailable
  if isempty(tmpD),
    tmpD=repmat(0,[mani AP.nLFPCh]);
  end
  for i=1:AP.nLFPCh
    eval(['tmpD(:,i)=strmread([AP.strmDir ''\'' rawCh(AP.LFPInd(i)).' strmType 'Fn],''nPts'',mani,''verbose'',0);']);
    if strcmpi(strmType,'gammaEnv') || strcmpi(strmType,'gammaNarrowEnv')
      % ** filter at theta frequency
      tmpD(:,i)=bafi(tmpD(:,i),WP.osi,AP.thetaCFreq,'rs',AP.rs);
    end
  end
  clear tmpd;
end

for i=1:length(r)
  if ~isempty(r(i).iPts)
    % try to find file containing user-defined attractor lags (defExpLag) and
    % check its contents
    isDefinedExpLag=false;
    defExpLag=[];
    lagFn=[DS.abfFn '_userdef_' AP.segmentType{i,1} '_' strmType 'cclag.mat'];
    if exist(lagFn,'file')
      load(lagFn);
      if ~exist('defExpLag','var')
        warndlg([lagFn ' must contain variable ''defExpLag'' for user-defined peak assignment - switching to automatic peak detection']);
      else
        [n1,n2]=size(defExpLag);
        if n1==AP.nAllLFPCh && n2==AP.nAllLFPCh
          % possibly implement other checks here in the future
          isDefinedExpLag=true;
          % ** cut down to expLag as used in code below, holding only values
          % from LFP channels TO BE ANALYZED
          defExpLag=defExpLag(AP.LFPccInd,AP.LFPccInd);
        else
          warndlg(['lag matrix as defined in ' lagFn ' should have dimension [' ...
            int2str(AP.nAllLFPCh)*[1 1] '] but in fact has dimension [' int2str([n1 n2]) ']']);
          defExpLag=[];
        end
      end
    end
    if isDefinedExpLag
      disp(['** user-defined lag values will be used for peak detection']);
    else
      disp(['** automatic peak detection']);
    end
    % assign preallocated templates to generic variables, omitting missing channels
    % (don't move two lines below, variables tmpccxxxx will be assigned full sized
    % templates again towards the end of current if case; doing this within loop saves memory)
    tmpccTemplate=WP.ccTemplate(AP.LFPccInd,AP.LFPccInd);
    tmpccDerTemplate=WP.ccDerTemplate(AP.LFPccInd,AP.LFPccInd);
    cccMn=tmpccTemplate;
    cccStd=tmpccTemplate;
    cccPeak=tmpccDerTemplate;          
    cccPeakMn=tmpccDerTemplate;          
    cccPeakStd=tmpccDerTemplate;          
    cccPeakT=tmpccDerTemplate;          
    cccPeakTMn=tmpccDerTemplate;          
    cccPeakTStd=tmpccDerTemplate;
    if strcmpi(strmType,'theta')
      cccPeakDecay=tmpccDerTemplate;
      cccPeakDecayMn=tmpccDerTemplate;
      cccPeakDecayStd=tmpccDerTemplate;
      cccPosPeakDecay=tmpccDerTemplate;
      cccPosPeakDecayMn=tmpccDerTemplate;
      cccPosPeakDecayStd=tmpccDerTemplate;
      cccNegPeakDecay=tmpccDerTemplate;
      cccNegPeakDecayMn=tmpccDerTemplate;
      cccNegPeakDecayStd=tmpccDerTemplate;
    end
    % intermediate results variables:
    % 2D matrix holding expected/preset lags between channels
    if isDefinedExpLag
      expLag=defExpLag;
      expAmp=[];
    else
      expLag=zeros(AP.nLFPCh);
      % 2D matrix holding corresponding amplitudes of chosen peaks (needed to weigh
      % computed lags)
      expAmp=expLag;
    end
    % segment-wise cc
    segCC=repmat(nan,[2*AP.ccLagPts+1 r(i).ni]);
    % segment-wise single parameter (peak decay)
    ninanny=repmat(nan,r(i).ni,1);
    % *** loop over diagonals ***
    for ki=1:nDiags
      ccix=WP.diagNccix{ki};
      nPairs=size(ccix,1);
      % set up the variable holding all peak amplitudes and lags of channel pairs
      % of current diagonal: as many columns as segments, as many rows as channel pairs;
      % each cell holding a 2 column array; 1st col peak times, 2nd col peak
      % amplitudes (note that this deletes former contents of tmpPeaks, if it
      % existed already)
      tmpPeaks=cell(nPairs,r(i).ni);

      % I. 
      % 1st loop over pairs in diagonals: 
      %    - computation of raw CC 
      %    - retrieval of ALL peak amplitudes and lags
      %    - estimation of expLag for current diagonal based on values in previous diagonal 
      for j=1:nPairs
        ccix1=ccix(j,1);
        ccix2=ccix(j,2);
        ci1=AP.LFPccInd(ccix1);
        ci2=AP.LFPccInd(ccix2);
        disp([r(i).segmentType ': XC ' strmType ' ' rawCh(AP.LFPInd(ccix1)).nm ' vs ' rawCh(AP.LFPInd(ccix2)).nm ': computing CC']);
        % load data ?
        if isempty(tmpD)
          % 1st channel
          tmpd1=eval(['strmread([AP.strmDir ''\'' rawCh(AP.LFPInd(ccix1)).' strmType 'Fn],''nPts'',r(i).iPts(end,2),''verbose'',0); ']);
          % load data (2nd channel)
          tmpd2=eval(['strmread([AP.strmDir ''\'' rawCh(AP.LFPInd(ccix2)).' strmType 'Fn],''nPts'',r(i).iPts(end,2),''verbose'',0); ']);
          if strcmpi(strmType,'gammaEnv') || strcmpi(strmType,'gammaNarrowEnv')            
            tmpd1=bafi(tmpd1,WP.osi,AP.thetaCFreq,'rs',AP.rs);
            tmpd2=bafi(tmpd2,WP.osi,AP.thetaCFreq,'rs',AP.rs);
          end
        else
          tmpd1=tmpD(:,ccix1);
          tmpd2=tmpD(:,ccix2);
        end
        if ki==1
          % * main diagonal: autocorr *
          % loop over intervals
          for k=r(i).ni:-1:1
            tmpIdx=r(i).iPts(k,1):r(i).iPts(k,2);
            segCC(:,k)=xxcorr(detrend(tmpd1(tmpIdx),'constant'),AP.ccLagPts,AP.ccScaleOpt);
          end
          % expLag is fine (0 along main diagonal)
        else 
          for k=r(i).ni:-1:1
            tmpIdx=r(i).iPts(k,1):r(i).iPts(k,2);
            segCC(:,k)=xxcorr(detrend(tmpd1(tmpIdx),'constant'),detrend(tmpd2(tmpIdx),'constant'),AP.ccLagPts,AP.ccScaleOpt);
          end
          % all but main diagonal: compute expected lag if not predefined
          if ki==2 
            if ~isDefinedExpLag
              % neighboring left element in same row + element in same column and diagonal of order 1 (ki==2))
              expLag(ccix1,ccix2)=expLag(ccix1,ccix2-1)+expLag(ccix2-1,ccix2);
            end
          else
            if ~isDefinedExpLag
              %             % below is the outline of an algorithm which may be useful because
              %             it integrates lag information from all channel combinations
              %             Mar 12 2005
              %             submatrix, side length: ki
              %             subMat=expLag(ccix1:ccix2,ccix1:ccix2);
              %             deLag=[];
              %             for mdi=ki-2:-1:1
              %               deLag(mdi)=sum(diag(subMat,mdi))-sum(diag(subMat,mdi-1));
              %             end
              %             expLag(ccix1,ccix2)=mean(deLag);
              wh=expAmp(ccix1,ccix1+1:ccix2-1)'+expAmp(ccix1+1:ccix2-1,ccix2);
              wh=(ki-2)*wh/sum(wh);
              expLag(ccix1,ccix2)=mean((expLag(ccix1,ccix1+1:ccix2-1)'+expLag(ccix1+1:ccix2-1,ccix2)).*wh);
            end
          end
        end % if:main diagonal
        % averaged CC  **segCC is not emptied**
        cccMn{ccix1,ccix2}(1:end,1)=mean(segCC,2);
        cccStd{ccix1,ccix2}(1:end,1)=std(segCC,0,2);
        % find all peaks 
        tmpr=evdeal(segCC(cci,:),WP.osi,'allpeaks');
        % alas, we need a loop to 
        % (i)  assemble peak times with amplitudes (**peak times are corrected
        %      for offset**)
        % (ii) determine peak decay (theta only)
        if strcmpi(strmType,'theta')
          % decay of AC as defined by 
          % a) slope of linear fit through log(abs(peaks))
          cccPeakDecay{ccix1,ccix2}=ninanny;
          % b) ratio of second largest and largest peaks
          cccPosPeakDecay{ccix1,ccix2}=ninanny;
          cccNegPeakDecay{ccix1,ccix2}=ninanny;
          for cpi=1:r(i).ni
            % tmpPeaks is not needed in this paragraph (decay times) but in
            % the following one
            tmpPeaks{j,cpi}=cat(2,tmpr.posPeakT{cpi}+tOffs,tmpr.posPeak{cpi});
            % variant a):
            % - combine pos and neg peaks
            tmpPNPeak=cat(1,tmpPeaks{j,cpi},cat(2,tmpr.negPeakT{cpi}+tOffs,tmpr.negPeak{cpi}));
            % - purge all with negative lags
            tmpPNPeak=tmpPNPeak(tmpPNPeak(:,1)>=0,:);
            if length(tmpPNPeak)>=2 & all(isfinite(tmpPNPeak))
              % - abs & log
              tmpPNPeak(:,2)=log(abs(tmpPNPeak(:,2)));
              % - fit and keep slope only 
              tmpb=regress(tmpPNPeak(:,2),[ones(size(tmpPNPeak,1),1) tmpPNPeak(:,1)]);
              cccPeakDecay{ccix1,ccix2}(cpi)=tmpb(2);
            end
            % variant b): sort according to amplitude and divide second
            % largest by largest
            tmpSortP=sort(tmpr.posPeak{cpi});
            if length(tmpSortP)>=2 & all(isfinite(tmpSortP))
              cccPosPeakDecay{ccix1,ccix2}(cpi)=tmpSortP(end-1)/tmpSortP(end);
            end
            tmpSortP=sort(tmpr.negPeak{cpi});
            if length(tmpSortP)>=2 & all(isfinite(tmpSortP))
              cccNegPeakDecay{ccix1,ccix2}(cpi)=tmpSortP(2)/tmpSortP(1);
            end
          end
          OKix=isfinite(cccPeakDecay{ccix1,ccix2});
          cccPeakDecayMn{ccix1,ccix2}=mean(cccPeakDecay{ccix1,ccix2}(OKix));
          cccPeakDecayStd{ccix1,ccix2}=std(cccPeakDecay{ccix1,ccix2}(OKix));
          OKix=isfinite(cccPosPeakDecay{ccix1,ccix2});
          cccPosPeakDecayMn{ccix1,ccix2}=mean(cccPosPeakDecay{ccix1,ccix2}(OKix));
          cccPosPeakDecayStd{ccix1,ccix2}=std(cccPosPeakDecay{ccix1,ccix2}(OKix));
          OKix=isfinite(cccNegPeakDecay{ccix1,ccix2});
          cccNegPeakDecayMn{ccix1,ccix2}=mean(cccNegPeakDecay{ccix1,ccix2}(OKix));
          cccNegPeakDecayStd{ccix1,ccix2}=std(cccNegPeakDecay{ccix1,ccix2}(OKix));
        else
          for cpi=1:r(i).ni
            tmpPeaks{j,cpi}=cat(2,tmpr.posPeakT{cpi}+tOffs,tmpr.posPeak{cpi});
          end
        end %strcmpi(strmType,'theta')
      end % for:j=1:nPairs
      
      % II. 
      % - search for the most likely 'path' along the mountain range
      %   defined by set of nth-neighbor CC
      % - update of expLag for current diagonal
      if ki==1
        % autocorr: path is along center peaks
      elseif ~isDefinedExpLag
        expLagDiag=diag(expLag,ki-1);
        % extract diagonal of matrix of averaged CC and cut down
        % tmpdd=diag(cccMn,ki-1); %original
        tmpdd = celldiag(cccMn, ki-1);
        tmpdd=cat(2,tmpdd{:});
        tmpdd=tmpdd(cci,:);
        % determine all peaks and, for the very very rare cases of no peak, max 
        tmpr=evdeal(tmpdd,WP.osi,{'allpeaks','minmax'});
        pix=cell(1,nPairs);
        nPaths=1;
        maxNPeak=1;
        spLag=[];
        spAmp=[];
        for j=1:nPairs
          % --- extract peaks:
          % 0. no peak found? This deplorable, reflecting too short a CC interval
          % length; however, don't halt program but put out max
          if isnan(tmpr.posPeak{j})
            warning('>>> no peak found: CC interval may be too short; substituting by max');
            tmpr.posPeak{j}=tmpr.max(j);
            tmpr.posPeakT{j}=tmpr.maxT(j);          
          end
          % 1. the number to be extracted depends on the amplitude of the max peak: the lower
          % it is, the more likely the current electrode pair spans a large distance and/or 
          % there may be a phase jump. In these cases it is important to retain many
          % peaks because either of the many small ones may be the 'right' one (based on
          % criteria defined below). On the other hand, if the max peak is large things are
          % almost always pretty clear, so let's just keep the max and a few of the largest 
          % side peaks
          [tmpP,tmpix]=sort(tmpr.posPeak{j}/max(tmpr.posPeak{j}),1);
          tmpP=flipud(tmpP);
          tmpix=flipud(tmpix);
          if max(tmpr.posPeak{j})>.75
            % maximally 2, all > .5*max peak amplitude
            tmpix=tmpix(tmpP>=.5);
            pix{j}=tmpix(1:min(2,length(tmpix)));
          elseif max(tmpr.posPeak{j})>.5
            % maximally 3, all > .2*max peak amplitude
            tmpix=tmpix(tmpP>=.2);
            pix{j}=tmpix(1:min(3,length(tmpix)));
          else
            % maximally 4, no limits on max peak amplitude
            pix{j}=tmpix(1:min(4,length(tmpix)));
          end
          maxNPeak=max(maxNPeak,length(pix{j}));
          nPaths=nPaths*length(pix{j});
          % as long as nPaths==1 keep collecting the actual lags so that we
          % have the values already in case there should be just one path
          if nPaths==1
            spLag(j,1)=tmpr.posPeakT{j}(pix{j})+tOffs;
            spAmp(j,1)=tmpr.posPeak{j}(pix{j});            
          end
        end
        if nPaths==1
          disp(['only one path found']);
          % update expLag & expAmp
          expLag=expLag+diag(spLag-diag(expLag,ki-1),ki-1);
          expLagDiag=diag(expLag,ki-1);
          expAmp=expAmp+diag(spAmp-diag(expAmp,ki-1),ki-1);          
        else
          spLag=[];          
          spAmp=[];                    
          disp(['choosing among ' int2str(nPaths) ' paths..']);
          % --- realizations of paths:
          % paths is a 3D matrix: the different paths in columns,
          % 1st slice the lags, second slice amplitudes
          paths=repmat(nan,[nPairs nPaths 2]);
          % unique lags and amplitudes (possibly needed for plot below)
          uLag=repmat(nan,maxNPeak,nPairs);
          uAmp=uLag;
          currLag=[];
          % a measure of how well expLag predicted the lags found (D=distance to expLag): 
          % [D(nearest peak)]/[D(nearest peak) + D(nearest peak on other side of expLag)]
          % values close to 0 mean good agreement between predicted and detected lag
          fracd=zeros(nPairs,2);
          nBlockRep=1;
          for j=1:nPairs
            nBlockRep=nBlockRep*length(pix{j});
            blockLen=nPaths/nBlockRep;
            % all peak lags found for the current pair
            currLag=tmpr.posPeakT{j}(pix{j})+tOffs;
            % peak lags
            paths(j,:,1)=reshape(repmat((currLag)',blockLen,nBlockRep/length(pix{j})),1,nPaths);
            % peak amplitudes
            paths(j,:,2)=reshape(repmat((tmpr.posPeak{j}(pix{j}))',blockLen,nBlockRep/length(pix{j})),1,nPaths);
            % unique lags and peaks
            uLag(1:length(currLag),j)=currLag;
            uAmp(1:length(currLag),j)=tmpr.posPeak{j}(pix{j});
            % extract closest peaks on either side of expLag, if any
            if any(currLag>=expLagDiag(j))
              fracd(j,1)=min(currLag(currLag>=expLagDiag(j))-expLagDiag(j));
            end
            if any(currLag<expLagDiag(j))            
               fracd(j,2)=min(expLagDiag(j)-currLag(currLag<expLagDiag(j)));
            end
          end
          % sort such that min distance is in column 1
          fracd=sort(fracd,2);
          % prevent the 'divide by zero' warning due to direct hits on single peak by setting second 
          % column arbitrarily to 1:
          fracd(~any(fracd,2),2)=1;
          % fractional distance
          fracd=fracd(:,1)./sum(fracd,2);

          % ******* now determine optimal path *******
          % 0. all criteria will always be computed even if they're not needed for the decision on the  
          %    optimal path because they will appear in plots
          % 1. check & possibly adjust current expLag diagonal: if prediction of lags by
          %    expLag is really bad (fracd >.4) do not rely on the criterion quantifying
          %    deviations from the expected path (crit1 below) but regard solely the minimal
          %    sum of squared lag jumps (crit3)
          
          % - crit3: sum of squared lag jumps (independent of expLag): 
          crit3=sum((diff(paths(:,:,1),1,1)).^2,1);
          % normalize to 1
          mc3=max(crit3);
          if mc3
            crit3=crit3/(mc3);
          end
          % check for plausibility of expLag:       
          badPix=find(fracd>maxFracd);
          if ~isempty(badPix)
            % should this happen with nearest neighbors rely on amplitude criterion
            if ki==2
              disp(['>>> bad prediction of pairs #' int2str(makerow(badPix)) ' by expLag; adjusting according to max average path altitude']);
              % .. will be done further below
            else
              goodPix=setdiff(1:nPairs,badPix);
              if isempty(goodPix)
                if nPairs==1
                  disp('>>> bad prediction of last-order neighbor lag by expLag; adjustment not possible');
                else
                  disp('>>> bad prediction of all lags by expLag; adjustment not possible');                  
                end
              else
                disp(['>>> bad prediction of pairs #' int2str(makerow(badPix)) ' by expLag; adjusting according to minimal jumpiness']);
                % find subpath (entity of all well-predicted lags) with least deviation 
                % from the subpath defined by lags in expLag
                crit0=sum(abs(paths(goodPix,:,1)-repmat(expLagDiag(goodPix),1,nPaths)),1);
                [pl,optix]=min(crit0);
                crit0=[];
                % now identify all full paths with this optimal subpath
                haix=find(~any(paths(goodPix,:,1)-repmat(paths(goodPix,optix,1),1,nPaths),1));
                % among those, determine the one with the least jitter
                [pl,optix]=min(crit3(haix));
                % don't forget..
                optix=haix(optix);
                % update expLag and expAmp
                expLag=expLag+diag(paths(:,optix,1)-diag(expLag,ki-1),ki-1);
                expLagDiag=diag(expLag,ki-1);
                expAmp=expAmp+diag(paths(:,optix,2)-diag(expAmp,ki-1),ki-1);
              end
            end
          end
          % now the other criteria: 
          % - little deviation from the path defined by lags in expLag
          crit1=sum(abs(paths(:,:,1)-repmat(diag(expLag,ki-1),1,nPaths)),1);
          % normalize to 1 (for plotting purposes): path with largest deviation is assigned value 1
          crit1=crit1/(max(crit1));
          % - high value of summed amplitudes
          crit2=sum(paths(:,:,2),1);
          % normalize to 1: divide by path with largest summed amplitude
          crit2=crit2/(max(crit2));
          if ki==2
            % first diagonal above main: this is the diagonal of nearest neighbor lags
            % which initializes prediction of lags between higher-order neighbors. 
            % Evaluate the optimal path based on amplitude criterion because we want it 
            % to wind along the crest which starts off with the highest amplitudes.
            [pl,optix]=max(crit2);
          else
            % if there is a choice among 10 or more paths, reject those whose normalized 
            % sum of squared lag jumps is way worse than that of the current expLag
            if nPaths>=10;
              % way worse: normalized sum > third of the way between expLag and 1 (the worst)
              expLagCrit3=sum(abs(diff(diag(expLag,ki-1))))/mc3;
              haix=find(crit3<=(1+expLagCrit3)/3);
              disp(['rejected ' int2str(nPaths-length(haix)) ' jittery paths']);
              % sometimes, things really get messed up and not a single path
              % qualifies
              if isempty(haix)
                haix=1:nPaths;
                disp(['rejected zero jittery paths because of very bad prediction by expLag']);                
              end
            else
              haix=1:nPaths;
            end
            % among those left, pick the one with the least deviation from expected path
            [pl,optix]=sort(crit1(haix),2);
            % ** in some cases more than one path has a minimum deviation from expected
            % path - make sure among those the one with the smallest sum of squared lag jumps
            % (crit3) is chosen **
            optix=optix(pl==pl(1));
            [pl,soptix]=min(crit3(haix(optix)));
            % don't forget..
            optix=haix(optix(soptix));
          end
          % update expLag and expAmp with values of currently found optimal path
          expLag=expLag+diag(paths(:,optix,1)-diag(expLag,ki-1),ki-1);
          expAmp=expAmp+diag(paths(:,optix,2)-diag(expAmp,ki-1),ki-1);          
          % plotting business
          ftitl=['order-' int2str(ki-1) ' neighbor CC'];
          figure(fh_ccLag);
          if nPaths<2
            clf
            text(.4,.4,'one path only');
            title(ftitl);
          else
            % anywhere between 500 and 1000 paths on plot (+ the optimal one)
            ppix=[optix 1:floor(max(ceil(nPaths/1000)*1000,1000)/1000):nPaths];
            if length(ppix)<nPaths+1, 
              titl='subset of paths';
            else
              titl='all paths';
            end
              
            subplot(2,4,2), 
            h1=plot(crit1(ppix)',crit2(ppix)','b.'); hold on
            h2=plot(crit1(optix),crit2(optix),'ro');            
            nicexyax; hold off;
            xlabel('lag deviation');
            ylabel('height');
            title(titl);
            
            subplot(2,4,3), 
            h1=plot(crit1(ppix)',crit3(ppix)','b.'); hold on
            h2=plot(crit1(optix),crit3(optix),'ro');            
            nicexyax; hold off;
            xlabel('lag deviation');
            ylabel('lag jitter');
            title(titl);

            subplot(2,4,4), 
            h1=plot(crit2(ppix)',crit3(ppix)','b.'); hold on
            h2=plot(crit2(optix),crit3(optix),'ro');            
            nicexyax; hold off;
            xlabel('height');
            ylabel('lag jitter');
            title(titl);

            subplot(2,4,6), 
            bar(fracd);
            lh=line(get(gca,'xlim'),[1 1]*maxFracd);
            set(lh,'color','m');
            set(gca,'ylim',[0 .5]);
            ylabel('nrmlzd pairwise lag dev');

            subplot(1,4,1), cla
            plotOffs=(1:nPairs)*max(diag(expAmp,ki-1))*1.2;
            set(gca,'color',get(gcf,'color'));
            hold on
            plot((discrete2cont(ccw(1):ccw(2),WP.osi*.001))',...
                 tmpdd+ones(size(tmpdd,1),1)*(plotOffs),'b-');
            plot(uLag,uAmp+repmat(plotOffs,maxNPeak,1),'ko');
            h=plot(paths(:,optix,1),paths(:,optix,2)+plotOffs','r-o');
            set(h,'linewidth',1.0');            
            niceyax;
            xlabel('lag (ms)');
            title(ftitl);
            hold off
          end
          drawnow;
        end % if:more than one path found
      end % if:main diagonal (autocorr)

      % III.
      % - for each segment, from all peaks previously found pull out the one closest to
      %   the corresponding expLag
      for j=1:nPairs
        ccix1=ccix(j,1);
        ccix2=ccix(j,2);
        ci1=AP.LFPccInd(ccix1);
        ci2=AP.LFPccInd(ccix2);
        disp([r(i).segmentType ': XC ' strmType ' ' rawCh(AP.LFPInd(ccix1)).nm ' vs ' rawCh(AP.LFPInd(ccix2)).nm ': pulling peaks']);
        tmppa=[];
        tmppl=[];
        % autocorr: pick max peak
        if ki==1        
          for cpi=r(i).ni:-1:1
            tmppa(cpi)=max(tmpPeaks{j,cpi}(:,2));
          end
          % peaks: amplitudes, mean & std
          cccPeak{ccix1,ccix2}=tmppa;
          cccPeakMn{ccix1,ccix2}=mean(cccPeak{ccix1,ccix2});
          cccPeakStd{ccix1,ccix2}=std(cccPeak{ccix1,ccix2});
          % peaks: latencies, mean & std
          cccPeakT{ccix1,ccix2}=zeros(size(tmppa));
          cccPeakTMn{ccix1,ccix2}=0;
          cccPeakTStd{ccix1,ccix2}=0;
        else
          for cpi=r(i).ni:-1:1
            [garnix,idefix]=min(abs(tmpPeaks{j,cpi}(:,1)-expLag(ccix1,ccix2)));
            tmppl(cpi)=tmpPeaks{j,cpi}(idefix,1);
            tmppa(cpi)=tmpPeaks{j,cpi}(idefix,2);
          end
          OKix=isfinite(tmppa);
          % peaks: amplitudes, mean & std (include nans)
          cccPeak{ccix1,ccix2}=tmppa;
          cccPeakMn{ccix1,ccix2}=mean(tmppa(OKix));
          cccPeakStd{ccix1,ccix2}=std(tmppa(OKix));
          % peaks: latencies, mean & std - offset is already subtracted
          cccPeakT{ccix1,ccix2}=tmppl;
          cccPeakTMn{ccix1,ccix2}=mean(tmppl(OKix));
          cccPeakTStd{ccix1,ccix2}=std(tmppl(OKix));
        end
      end
    end % for: # diagonals

    % embed all reduced generic 2D CC variables in templates with original dimension
    tmpccTemplate=WP.ccTemplate;
    tmpccTemplate(AP.LFPccInd,AP.LFPccInd)=cccMn;
    eval(['r(i).' STshort 'CCMn=tmpccTemplate;']);    
    
    tmpccTemplate=WP.ccTemplate;
    tmpccTemplate(AP.LFPccInd,AP.LFPccInd)=cccStd;
    eval(['r(i).' STshort 'CCStd=tmpccTemplate;']);    
    
    tmpccDerTemplate=WP.ccDerTemplate;    
    tmpccDerTemplate(AP.LFPccInd,AP.LFPccInd)=cccPeak;    
    eval(['r(i).' STshort 'CCPeak=tmpccDerTemplate;']);    

    tmpccDerTemplate=WP.ccDerTemplate;    
    tmpccDerTemplate(AP.LFPccInd,AP.LFPccInd)=cccPeakMn;    
    eval(['r(i).' STshort 'CCPeakMn=tmpccDerTemplate;']);    

    tmpccDerTemplate=WP.ccDerTemplate;    
    tmpccDerTemplate(AP.LFPccInd,AP.LFPccInd)=cccPeakStd;    
    eval(['r(i).' STshort 'CCPeakStd=tmpccDerTemplate;']);    

    tmpccDerTemplate=WP.ccDerTemplate;    
    tmpccDerTemplate(AP.LFPccInd,AP.LFPccInd)=cccPeakT;    
    eval(['r(i).' STshort 'CCPeakT=tmpccDerTemplate;']);    

    tmpccDerTemplate=WP.ccDerTemplate;    
    tmpccDerTemplate(AP.LFPccInd,AP.LFPccInd)=cccPeakTMn;    
    eval(['r(i).' STshort 'CCPeakTMn=tmpccDerTemplate;']);    

    tmpccDerTemplate=WP.ccDerTemplate;    
    tmpccDerTemplate(AP.LFPccInd,AP.LFPccInd)=cccPeakTStd;    
    eval(['r(i).' STshort 'CCPeakTStd=tmpccDerTemplate;']);    
    
    if strcmpi(strmType,'theta')

      tmpccDerTemplate=WP.ccDerTemplate;
      tmpccDerTemplate(AP.LFPccInd,AP.LFPccInd)=cccPeakDecay;
      eval(['r(i).' STshort 'CCPeakDecay=tmpccDerTemplate;']);

      tmpccDerTemplate=WP.ccDerTemplate;
      tmpccDerTemplate(AP.LFPccInd,AP.LFPccInd)=cccPeakDecayMn;
      eval(['r(i).' STshort 'CCPeakDecayMn=tmpccDerTemplate;']);

      tmpccDerTemplate=WP.ccDerTemplate;
      tmpccDerTemplate(AP.LFPccInd,AP.LFPccInd)=cccPeakDecayStd;
      eval(['r(i).' STshort 'CCPeakDecayStd=tmpccDerTemplate;']);

      tmpccDerTemplate=WP.ccDerTemplate;
      tmpccDerTemplate(AP.LFPccInd,AP.LFPccInd)=cccPosPeakDecay;
      eval(['r(i).' STshort 'CCPosPeakDecay=tmpccDerTemplate;']);

      tmpccDerTemplate=WP.ccDerTemplate;
      tmpccDerTemplate(AP.LFPccInd,AP.LFPccInd)=cccPosPeakDecayMn;
      eval(['r(i).' STshort 'CCPosPeakDecayMn=tmpccDerTemplate;']);

      tmpccDerTemplate=WP.ccDerTemplate;
      tmpccDerTemplate(AP.LFPccInd,AP.LFPccInd)=cccPosPeakDecayStd;
      eval(['r(i).' STshort 'CCPosPeakDecayStd=tmpccDerTemplate;']);
      
      tmpccDerTemplate=WP.ccDerTemplate;
      tmpccDerTemplate(AP.LFPccInd,AP.LFPccInd)=cccNegPeakDecay;
      eval(['r(i).' STshort 'CCNegPeakDecay=tmpccDerTemplate;']);

      tmpccDerTemplate=WP.ccDerTemplate;
      tmpccDerTemplate(AP.LFPccInd,AP.LFPccInd)=cccNegPeakDecayMn;
      eval(['r(i).' STshort 'CCNegPeakDecayMn=tmpccDerTemplate;']);

      tmpccDerTemplate=WP.ccDerTemplate;
      tmpccDerTemplate(AP.LFPccInd,AP.LFPccInd)=cccNegPeakDecayStd;
      eval(['r(i).' STshort 'CCNegPeakDecayStd=tmpccDerTemplate;']);
    end
  end % if ~isempty(r(i).iPts)
end % for i=1:length(r)
% clean up
if ishandle(fh_ccLag), delete(fh_ccLag); end
% reset graphics settings
labelscale('fontSz',7,'scaleFac',1.0,'lineW',1.0,'markSz',4); 

warning('on','stats:regress:RankDefDesignMat');