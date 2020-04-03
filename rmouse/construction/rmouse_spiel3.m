pack
switch strmType
  case 'delta'
    % detect peaks within this interval centered around t=0
    ccw=cont2discrete(AP.deccw,osi*.001,'intv',0)-1; % sic!
  case 'theta'
    ccw=cont2discrete(AP.thccw,osi*.001,'intv',0)-1; % sic!
  case 'thetaHiEnv'
    ccw=cont2discrete(AP.deccw,osi*.001,'intv',0)-1; % sic!
  case 'thetaLoEnv'
    ccw=cont2discrete(AP.ccLagPts*[-.95 .95],osi*.001,'intv',0)-1; % sic!
  case 'gamma'
    ccw=cont2discrete(AP.gaccw,osi*.001,'intv',0)-1; % sic!
  case 'gammaEnv'
    ccw=cont2discrete(AP.gaccw,osi*.001,'intv',0)-1; % sic!
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
tOffs=discrete2cont(ccw(1),osi*.001,'intv',0);
% the number of diagonals in cc matrix
nDiags=length(diagNccix);

% --- load it up, Scotty (if you can)
% if RAM permits, upload all channels completely here (thus circumventing multiple 
% read operations with each segment type and channel)
mame=max(cat(1,r(:).dmem));
mani=max(cat(1,r(:).lastPt));
if mame*(nLFPCh+2)/nLFPCh<=WP.maxRAM;
  tmpD=repmat(0,[mani nLFPCh]);
  for i=1:nLFPCh
    eval(['tmpD(:,i)=strmread([AP.strmDir ''\'' rawCh(AP.LFPInd(i)).' strmType 'Fn],''nPts'',mani,''verbose'',0);']);
    if strcmpi(strmType,'gammaEnv')
      % ** filter at theta frequency
      tmpD(:,i)=bafi(tmpD(:,i),osi,AP.thetaCFreq,'rs',AP.rs);
    end
  end
  clear tmpd;
else 
  tmpD=[];
end

for i=1:length(r)
  if ~isempty(r(i).iPts)
    % assign preallocated templates to generic variables, omitting missing channels
    % (don't move two lines below, variables tmpccxxxx will be assigned full sized
    % templates again towards the end of current if case; doing this within loop saves memory)
    tmpccTemplate=ccTemplate(AP.LFPccInd,AP.LFPccInd);
    tmpccDerTemplate=ccDerTemplate(AP.LFPccInd,AP.LFPccInd);
%     cccMn=tmpccTemplate;
%     cccStd=tmpccTemplate;    
%     cccPeak=tmpccDerTemplate;          
%     cccPeakMn=tmpccDerTemplate;          
%     cccPeakStd=tmpccDerTemplate;          
%     cccPeakT=tmpccDerTemplate;          
%     cccPeakTMn=tmpccDerTemplate;          
%     cccPeakTStd=tmpccDerTemplate;                    


    
    % intermediate results variables:
    % 2D matrix holding expected lags between channels
    expLag=zeros(nLFPCh);
    % 2D matrix holding corresponding amplitudes of chosen peaks (needed to weigh
    % computed lags)
    expAmp=expLag;
    % *** loop over segments ***
    for k=1:r(i).ni
      tmpIdx=r(i).iPts(k,1):r(i).iPts(k,2);
      % set up the variable holding all peak amplitudes and lags of channel pairs;
      % each cell holding a 2 column array; 1st col peak times, 2nd col peak
      % amplitudes (note that this deletes former contents of tmpPeaks, if it
      % existed already)
      tmpPeaks=cell(nLFPCh);
      % *** loop over diagonals ***
%      for ki=1:nDiags
      for ki=1:2
        ccix=diagNccix{ki};
        nPairs=size(ccix,1);
        ccc=repmat(nan,[2*AP.ccLagPts+1,nPairs]);
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
            eval([ 'tmpd1=strmread([AP.strmDir ''\'' rawCh(AP.LFPInd(ccix1)).' strmType 'Fn],''nPts'',r(i).iPts(end,2),''verbose'',0); ']);
            % load data (2nd channel)
            eval([ 'tmpd2=strmread([AP.strmDir ''\'' rawCh(AP.LFPInd(ccix2)).' strmType 'Fn],''nPts'',r(i).iPts(end,2),''verbose'',0); ']);
            if strcmpi(strmType,'gammaEnv')
              tmpd1=bafi(tmpd1,osi,AP.thetaCFreq,'rs',AP.rs);
              tmpd2=bafi(tmpd2,osi,AP.thetaCFreq,'rs',AP.rs);
            end
          else
            tmpd1=tmpD(:,ccix1);
            tmpd2=tmpD(:,ccix2);
          end
          if ki==1
            % * main diagonal: autocorr *
            ccc(:,j)=xxcorr(detrend(tmpd1(tmpIdx),'constant'),AP.ccLagPts,AP.ccScaleOpt);
            % expLag is fine (0 along main diagonal)
          else % if:main diagonal
            ccc(:,j)=xxcorr(detrend(tmpd1(tmpIdx),'constant'),detrend(tmpd2(tmpIdx),'constant'),AP.ccLagPts,AP.ccScaleOpt);
            % other diagonals: compute expected lag 
            if ki==2
              % neighboring left element in same row + element in same column and diagonal of order 1 (ki==2))
              expLag(ccix1,ccix2)=expLag(ccix1,ccix2-1)+expLag(ccix2-1,ccix2);
            else
              % weigh by CC peak amplitudes
              wh=expAmp(ccix1,ccix1+1:ccix2-1)'+expAmp(ccix1+1:ccix2-1,ccix2);
              wh=(ki-2)*wh/sum(wh);
              expLag(ccix1,ccix2)=mean((expLag(ccix1,ccix1+1:ccix2-1)'+expLag(ccix1+1:ccix2-1,ccix2)).*wh);
            end
          end % if:main diagonal
%           % find all peaks 
%           tmpr=evdeal(ccc(cci,j),osi,'allpeaks');
%           % alas, we need to assemble peak times with amplitudes
%           % **peak times are corrected for offset
%           tmpPeaks{ccix1,ccix2}=cat(2,tmpr.posPeakT{1}+tOffs,tmpr.posPeak{1});
        end % for:j=1:nPairs
      % II. 
      % - search for the most likely 'path' along the mountain range
      %   defined by set of nth-neighbor CC
      % - update of expLag for current diagonal
      if ki==1
        % autocorr: path is along center peaks
      else
        expLagDiag=diag(expLag,ki-1);
        tmpdd=ccc(cci,:);
        % determine all peaks and, for the very very rare cases of no peak, max 
        tmpr=evdeal(tmpdd,osi,{'allpeaks','minmax'});
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
            % maximally 5, no limits on max peak amplitude
            pix{j}=tmpix(1:min(5,length(tmpix)));
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
                % find subpath with least deviation from the subpath defined by lags in expLag
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
%          if nPaths<2
          if nPaths<1
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
            plot((discrete2cont(ccw(1)+1:ccw(2)+1,osi*.001)+osi*.001*.5)',...
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
    end
    end
   end

end
% clean up
if ishandle(fh_ccLag), delete(fh_ccLag); end
% reset graphics settings
labelscale('fontSz',7,'scaleFac',1.0,'lineW',1.0,'markSz',4); 
clear d ccc* tmp* npts i k chI* ptr nptr nix nixIdx mame mani
