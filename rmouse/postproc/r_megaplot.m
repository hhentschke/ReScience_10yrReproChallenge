function r_megaplot(mAP,varargin)
% ** function r_megaplot(mAP,varargin)
%    generates collection of plots from all results files for single variable.
%    Needs ANPAR and DSET (=concatenation of individual and matching(!) AP and DS)
%        >>> INPUT VARIABLES >>>
% NAME      TYPE/DEFAULT          DESCRIPTION
% mAP       char array            file name of the 'master' AP (analysis parameters) 
%                                 that had been used in computing the results to be
%                                 extracted and plotted by this function. May be set to empty
%                                 string ''  
% behav     cell array of chars   any of 'immobile', 'exploring'
% rv        cell array of chars   most of legal fields of r, like 'thCCMn' (a
%                                 waveform resulting in line or contour plot)
%                                 'thCCPeakMn' (scalars resulting in matrix
%                                 plots), etc.
% wavefm    scalar,0              by default, this script produces contour plots
%                                 or matrix plots, depending on the contents of
%                                 rv. If wavefm is nonzero, waveforms will be
%                                 plotted as lines instead
% pn        char array, []        if nonempty, this string will be inserted in file name
%                                 of plot (e.g. sth like 'wt' may be useful)
% figdir    char array, []        partial directory for figure files, e.g. '\beta3_wtko\figures'
% 

global ANPAR DSET AP DS

% ----- default values & varargin -----
behav={'exploring'};
rv={'thgaeCCMn'};
wavefm=0;
pn=[];
figdir=[];
printas={'-djpeg90'};{[]};
pvpmod(varargin);

behavType={...
    'exploring', [6     9],  [.5 1 .5],   'o';...
    'immobile',  [2.5   4],  [.6 .4 .1],  's';...
};
rvType={'deCCMn','thCCMn','gaCCMn','thLoeCCMn','thHieCCMn','gaeCCMn','detheCCMn','thgaeCCMn',...
  'deCCPeakMn','thCCPeakMn','gaCCPeakMn','thLoeCCPeakMn','thHieCCPeakMn','gaeCCPeakMn',...
  'deCCPeakStd','thCCPeakStd','gaCCPeakStd','thLoeCCPeakStd','thHieCCPeakStd','gaeCCPeakStd',...
  'deCCPeakTMn','thCCPeakTMn','gaCCPeakTMn','thLoeCCPeakTMn','thHieCCPeakTMn','gaeCCPeakTMn',...
  'deCCPeakTStd','thCCPeakTStd','gaCCPeakTStd','thLoeCCPeakTStd','thHieCCPeakTStd','gaeCCPeakTStd',...
  'fComod','fComodP','rawCohMn',...
  'rawDePEMn','rawThPEMn','rawThNarrowPEMn','rawBePEMn','rawGaPEMn','rawRiPEMn',...
  'rawgaeCohTh','rawgaeCohPeak','rawgaeCohPeakF',};

% --------- check of input vars
if exist('plotmd','var'),
  error('input variable ''plotmd'' does not exist anymore - see input var ''wavefm''');
end

nBehav=length(behav);
nRv=length(rv);


%       if isempty(strmatch(behav{bi},behavType(:,1)))
%         error('check input var ''behav''');
%       end
%       if isempty(strmatch(rv{rvi},rvType))
%         error('check input var ''rv''');
%       end


% --- check ANPAR & DSET + other preps
[n1 n2 n3]=size(ANPAR);
[dn1 dn2 dn3]=size(DSET);
if ~isequal([n1 n2 n3],[dn1 dn2 dn3]), error('DSET and ANPAR must have equal dimension'); end
if n1>1, error([mfilename ' does not deal with multi-row ANPARs and DSETs - split them up in separate data sets']); end
rmouse_ini;

% ---------- plot settings
close all;
labelscale('fontSz',6,'scaleFac',1.0,'lineW',1.0,'markSz',6);
% color associated with principal channel
pcCol=[1 .6 .6];
% # of levels in contour plots
nLev=20;
switch n1
  case 1
    if wavefm
      % depending on number of plots, design layout of subplots (approx 2:8 ratio)
      nRow=round(sqrt(2/8*n1*n2));
      nCol=ceil(n1*n2/nRow);
    else
      % depending on number of plots, design layout of subplots (approx 7:3 ratio)
      nRow=round(sqrt(7/3*n1*n2));
      nCol=ceil(n1*n2/nRow);
    end
  case 2
    error('plot layout undefined for n=2 rows of ANPAR');
  case 3
    error('plot layout undefined for n=3 rows of ANPAR');
  otherwise
    error('plot layout undefined for n>3 rows of ANPAR');
end


% ------------- the works
% loop over data sets:
% one experiment per column, concentration down the columns
for ci=1:n2
  for ri=1:n1
    AP=ANPAR(ri,ci);
    DS=DSET(ri,ci);
    % append/overwrite common settings, if any
    try
      if ~isempty(mAP)
        eval([mAP ';']);
      end
    catch
      lasterr
      warndlg(['There was a problem finding or reading the master AP (' mAP ') - if ' mfilename ' produces an error, this problem is the most likely cause of it']);
    end
    rmouse_apcheck;
    % --- channels
    rawCh=rmouse_chan;
    % some shorties
    pcIdx=AP.pcIdx;
    pcInd=AP.pcInd;
    LFPpcInd1=AP.LFPpcInd1;
    LFPpcInd2=AP.LFPpcInd2;
    LFPccInd=AP.LFPccInd;
    % template
    tempo=repmat(nan,nAllLFPCh,1);
    % --- data & paths
    % if dpath does not contain a drive letter, pre-pend WP.rootPath
    if isempty(strfind(DS.dpath,':')), DS.dpath=[WP.rootPath DS.dpath]; end
    % same with results path & stream dir
    if isempty(strfind(AP.resPath,':')), AP.resPath=[WP.rootPath AP.resPath]; end
    if exist([DS.dpath '\' DS.abfFn '.abf'],'file')
      [nix,nix2,abfi]=abfload([DS.dpath '\' DS.abfFn '.abf'],'info');      
      % abfi=abfinfo([DS.dpath '\' DS.abfFn '.abf'],'verbose',0);
      % sampling interval
      osi=abfi.si;
    else
      warndlg([DS.dpath '\' DS.abfFn ' does not exist; assuming sampling interval of 1000 us'])
    end
    % load results var..
    load([AP.resPath '\' AP.resFn],'r');
    for bi=1:nBehav
      % determine index to results obtained with specified behavior in r
      i=strmatch(behav{bi},{r(:).segmentType});
      pcol=behavType{strmatch(behav{bi},behavType(:,1)),3};
      psym=behavType{strmatch(behav{bi},behavType(:,1)),4};
      % subplot index
      switch n1
        case 1
          spix=(ci-1)*n1+ri;
        case 2
        case 3
        otherwise
      end
      % ** Instructions for transplanting code from rmouse_p to here **
      % --- delete:
      % for bi..
      % sph=subplot('position',[xpos(lct) ypos(row) xlen ylen]);
      % lct=lct+1;
      % pcol=..
      % cbh=colorbar...
      % title..
      % --- change
      % set(lh,'color',pcCol,'linestyle','--','linewidth',4);
      % --- add
      % colormap(bone);
      for rvi=1:nRv
        figName=[pn '_' behav{bi}(1:3) '_' rv{rvi}];
        mkfig(figName); orient tall
        subplot(nRow,nCol,spix);
        switch rv{rvi}
          case {'rawgaeCohPeak','rawgaeCohPeakF','rawgaeCohTh'}
            % colormap
            nColors=128;
            ccprop={'peak ampl','peak freq'};
            switch rv{rvi}
              case 'rawgaeCohPeak'
                cLim=[0 0.8];
                colormap(spring(nColors));
              case 'rawgaeCohPeakF'
                cLim=[7 10];
                colormap(summer(nColors));
              case 'rawgaeCohTh'
                cLim=[0 0.4];
                colormap(spring(nColors));
            end
            if isfield(r,[rv{rvi}])
              eval(['cohMat=r(i).' rv{rvi} ';']);
              if ~isempty(cohMat)
                cohMat=permute(cohMat,[3 2 1]);
                ih=imagesc(cohMat,cLim);
                cph=gca;
                set(cph,'xaxisloc','top',...
                  'ylim',[.5 nAllLFPCh+.5],'ytick',[1:nAllLFPCh],'yticklabel',WP.xtl2,...
                  'xlim',[.5 nAllLFPCh+.5],'xtick',[1:nAllLFPCh],'xticklabel',WP.xtl2);
                axis square
                % lines to indicate principal electrode
                lh=line(LFPpcInd2*[1 1]',[.5 nAllLFPCh+.5]);
                set(lh,'color',[.6 .6 .6],'linestyle','-');
                lh=line([.5 nAllLFPCh+.5],LFPpcInd2*[1 1]');
                set(lh,'color',[.6 .6 .6],'linestyle','-');
                xlabel('raw');
                ylabel('gammaEnv');
                if ci==n2 && ri==n1
                  cbh=colorbar;
                end
              end
            end
          
          % --------------------
          % --- the matrix plots
          % --------------------
          case {'deCCPeakMn','deCCPeakTMn','thCCPeakMn','thCCPeakTMn',...
              'gaCCPeakMn','gaCCPeakTMn','gaeCCPeakMn','gaeCCPeakTMn',...
              'thLoeCCPeakMn','thLoeCCPeakTMn','thHieCCPeakMn','thHieCCPeakTMn',...
              'deCCPeakStd','deCCPeakTStd','thCCPeakStd','thCCPeakTStd',...
              'gaCCPeakStd','gaCCPeakTStd','gaeCCPeakStd','gaeCCPeakTStd',...
              'thLoeCCPeakStd','thLoeCCPeakTStd','thHieCCPeakStd','thHieCCPeakTStd'};
            nColors=128;
            if strfind(rv{rvi},'PeakMn')
              colormap(jet(nColors));
            elseif strfind(rv{rvi},'PeakStd')
              colormap(gray(nColors));
            elseif strfind(rv{rvi},'PeakTMn')
              colormap(coma('amberturquois','ncols',nColors));
            elseif strfind(rv{rvi},'PeakTStd')
              colormap(gray(nColors));
            end
            ccMatTemplate=repmat(nan,[nAllLFPCh nAllLFPCh]);
            % set limits so that plot-distorting outliers can be kicked out 
            switch rv{rvi}
              case 'deCCPeakMn'
                maxV=1;
                minV=.0;
              case 'deCCPeakStd'
                maxV=.4;
                minV=.0;
              case 'deCCPeakTMn'
                maxV=150;
                minV=-150;
              case 'deCCPeakTStd'
                maxV=150;
                minV=0;

              case 'thCCPeakMn'
                maxV=1;
                minV=.0;
              case 'thCCPeakStd'
                maxV=.4;
                minV=.0;
              case 'thCCPeakTMn'
                maxV=75;
                minV=-30;
              case 'thCCPeakTStd'
                maxV=40;
                minV=0;

              case 'thLoeCCPeakMn'
                maxV=1;
                minV=.0;
              case 'thLoeCCPeakStd'
                maxV=.4;
                minV=.0;
              case 'thLoeCCPeakTMn'
                maxV=250;
                minV=-250;
              case 'thLoeCCPeakTStd'
                maxV=200;
                minV=0;

              case 'thHieCCPeakMn'
                maxV=1;
                minV=.0;
              case 'thHieCCPeakStd'
                maxV=.4;
                minV=.0;
              case 'thHieCCPeakTMn'
                maxV=250;
                minV=-250;
              case 'thHieCCPeakTStd'
                maxV=200;
                minV=0;

              case 'gaCCPeakMn'
                maxV=1;
                minV=.0;
              case 'gaCCPeakStd'
                maxV=.4;
                minV=.0;
              case 'gaCCPeakTMn'
                maxV=10;
                minV=-10;
              case 'gaCCPeakTStd'
                maxV=4;
                minV=0;

              case 'gaeCCPeakMn'
                maxV=1;
                minV=.0;
              case 'gaeCCPeakStd'
                maxV=.4;
                minV=.0;
              case 'gaeCCPeakTMn'
                maxV=10;
                minV=-10;
              case 'gaeCCPeakTStd'
                maxV=4;
                minV=0;
              otherwise
                error('illegal streamType in CC matrix plot');
            end

            % scale the data such that [min max] maps to [1 nColors]
            sFac=(nColors-1)/(maxV-minV);
            sOffs= -minV;
            if isfield(r,rv{rvi})
              eval(['cm1=diag(r(i).' rv{rvi} ');']);
              if length(cm1)>1
                eval(['cm1=r(i).' rv{rvi} ';']);
                ccMat=ccMatTemplate;
                % transfer values from cell array in array - has to be done elementwise because
                % there are empty cells
                for g=1:nAllLFPCh^2
                  if ~isempty(cm1{g})
                    ccMat(g)=cm1{g};
                  end
                end
                % deal with outliers (set to nan)
                ccMat((ccMat)>maxV)=nan;
                ccMat((ccMat)<minV)=nan;
                ccMat=round(sOffs+(sFac*ccMat))+1;
                % ---
                ii=1;
                ih=image(ccMat);
                set(ih,'CDataMapping','direct')
                cph=gca;
                set(cph,'xaxisloc','top',...
                  'ylim',[.5 nAllLFPCh+.5],'ytick',[1:nAllLFPCh],'yticklabel',WP.xtl2,...
                  'xlim',[.5 nAllLFPCh+.5],'xtick',[1:nAllLFPCh],'xticklabel',WP.xtl2);
                axis square
                % lines to indicate principal electrode
                lh=line(LFPpcInd2*[1 1]',[.5 LFPpcInd2]);
                set(lh,'color',[.6 .6 .6],'linestyle','-');
                lh=line([LFPpcInd2 nAllLFPCh+.5],LFPpcInd2*[1 1]');
                set(lh,'color',[.6 .6 .6],'linestyle','-');
                % inflate plot
                rexy('ax',cph,'xfac',1.1,'yfac',1.1);
              end
            end

            % -------------------------
            % --- comodulograms
            % -------------------------
          case 'fComod'
            colormap(flipud(coma('amberturquois','ncols',255)));
            % assume that bulk of correlation coeffs remain in that interval
            ccLim=[-.5 .5];
            if isfield(r,'fComod')
              cph=imagesc(r(1).comF,r(1).comF,r(i).fComod(:,:,AP.LFPpcInd1),ccLim);
              axis square
              set(gca,'xtick',[10:20:r(1).comF(end)],'ytick',[10:20:r(1).comF(end)],...
                'xaxislocation','top');
              grid on
            end

          case 'fComodP'
            colormap(gray);
            % assume that bulk of correlation coeffs remain in that interval
            pLim=[0 .05];
            if isfield(r,'fComod')
              cph=imagesc(r(1).comF,r(1).comF,r(i).fComodP(:,:,AP.LFPpcInd1),pLim);
              axis square
              set(gca,'xtick',[10:20:r(1).comF(end)],'ytick',[10:20:r(1).comF(end)],...
                'xaxislocation','top');
              grid on
            end
            % -------------------------
            % contour and/or line plots
            % -------------------------
          case 'rawCohMn'
            colormap(hot);
            % # of contours
            nLev=20;
            srfc=[];
            % large range (0-100 Hz)
            frix=find(r(1).F>=0 & r(1).F<=100);
            lct=1;
            hold on
            if length(diag(r(i).rawCohMn))>1
              % freq range from which average theta coherence had been calculated
              narrThF=r(i).rawPMnPeakT{AP.LFPpcInd2,AP.LFPpcInd2}+[-1 1];
              srfc=cat(3, srfc, repmat(nan,length(frix),nAllLFPCh));
              % the lines below take into account non-analyzed channels (which have
              % single nans in the corresponding cells)
              for g=1:nAllLFPCh
                ix1=min(g,AP.LFPpcInd2);
                ix2=max(g,AP.LFPpcInd2);
                if length(r(i).rawCohMn{ix1,ix2})>1
                  srfc(1:length(frix),g,lct)=r(i).rawCohMn{ix1,ix2}(frix);
                else
                  srfc(1:length(frix),g,lct)=nan;
                end
              end
              % --- contour plot
              [c,cph]=contourf(1:size(srfc,2),r(1).F(frix),srfc(:,:,lct),nLev);
              clear c;
              axis tight
              set(cph(:),'linestyle','none');
              set(gca,'clim',[0 1]);
              cph=gca;
              % matlab V. 7.0.1 completely fucks up overlaid axes with colorbars, so let it be
              % for now
              % cbh=colorbar(cbPos);
              % use the contour plot's x axis as electrode # indicator
              set(cph,'xaxisloc','top','yaxisloc','right','ylim',[r(1).F(frix(1)) r(1).F(frix(end))],'ytick',[],...
                'xlim',[1 nAllLFPCh],'xtick',[1:nAllLFPCh],'xticklabel',WP.xtl2);
              % x axis limits must be fixed for aligned overlay with surface plot
              ax2=axes('position',get(cph,'position'),'color','none');
              hold on;
              set(gca,'ylim',[r(1).F(frix(1)) r(1).F(frix(end))],'xlim',[WP.elx(1) WP.elx(end)],...
                'xtick',WP.elx,'xticklabel',WP.xtl1);
              axis manual;
              grid on;
              % line to indicate principal electrode
              lh=line(WP.elx(LFPpcInd2)*[1 1]',get(gca,'ylim'));
              set(lh,'color',pcCol,'linestyle',':');
              % push to background
              set(gca,'children',circshift(get(gca,'children'),-1))
              % lines to indicate freq ranges from which average coherence had been calculated
              lh=line([WP.elx(1) WP.elx(end)]'*[1 1],[1;1]*AP.delta,'color','c');
              lh=line([WP.elx(1) WP.elx(end)]'*[1 1],[1;1]*narrThF,'color','b');
              lh=line([WP.elx(1) WP.elx(end)]'*[1 1],[1;1]*AP.gamma,'color','g');
              title(['Coherence;' r(i).segmentType]);
              ylabel('Frequency (Hz)');
              grid on;
            end


          case 'thHieCCMn'
            if isfield(r,'thHieCCMn')
              colormap(bone);
              % thetaHi env (range to plot: the lesser of [+/- 1 slow theta period, AP.ccLagPts])
              ccw=[-1 1]*min(AP.ccLagPts,cont2discrete(1*1000/AP.theta(1),osi*.001,'intv',1));
              cci=AP.ccLagPts+[ccw(1):ccw(2)]+1;
              ccp=discrete2cont(cci-AP.ccLagPts,osi*.001,'intv',0);
              hold on
              if length(diag(r(i).thHieCCMn))>1
                srfc=repmat(nan,length(cci),nAllLFPCh);
                % the lines below take into account non-analyzed channels (which have
                % single nans in the corresponding cells)
                for g=1:nAllLFPCh
                  ix1=min(g,LFPpcInd2);
                  ix2=max(g,LFPpcInd2);
                  if length(r(i).thHieCCMn{ix1,ix2})>1
                    srfc(1:length(cci),g)=r(i).thHieCCMn{ix1,ix2}(cci);
                    % symmetry
                    if g>LFPpcInd2, srfc(:,g)=flipud(srfc(:,g)); end
                  else
                    srfc(1:length(cci),g)=nan;
                  end
                end
                % --- contour plot
                [c,cph]=contourf(srfc,nLev);
                clear c;
                axis tight
                set(cph(:),'linestyle','none');
                cph=gca;
                % use the contour plot's x axis as electrode # indicator
                set(cph,'xaxisloc','top','yaxisloc','right','ytick',[],...
                  'xlim',[1 nAllLFPCh],'xtick',[1:nAllLFPCh],'xticklabel',WP.xtl2);
                % x axis limits must be fixed for aligned overlay with surface plot
                ax2=axes('position',get(cph,'position'),'color','none');
                hold on;
                set(gca,'ylim',[ccp(1) ccp(end)],'xlim',[WP.elx(1) WP.elx(end)],...
                  'xtick',WP.elx,'xticklabel',WP.xtl1);
                axis manual;
                tmpMn=cat(1,r(i).thHieCCPeakTMn{1:LFPpcInd2,LFPpcInd2});
                % symmetry - don't forget to invert
                tmpMn=[tmpMn; -1*cat(1, r(i).thHieCCPeakTMn{LFPpcInd2,LFPpcInd2+1:end})];
                tmpStd=cat(1, r(i).thHieCCPeakTStd{1:LFPpcInd2,LFPpcInd2});
                tmpStd=[tmpStd; cat(1, r(i).thHieCCPeakTStd{LFPpcInd2,LFPpcInd2+1:end})];
                ph=errorbar(WP.elx,tmpMn,tmpStd,[psym '-']);
                set(ph,'color',pcol);
                grid on;
                % line to indicate principal electrode
                lh=line(WP.elx(LFPpcInd2)*[1 1]',get(gca,'ylim'));
                set(lh,'color',pcCol,'linestyle','--','linewidth',4);
                % push to background
                set(gca,'children',circshift(get(gca,'children'),-1))
              end
            end
          case 'thLoeCCMn'
            if isfield(r,'thLoeCCMn')
              colormap(bone);
              % thetaLo env (range to plot: the lesser of [+/- 1 slow theta period, AP.ccLagPts])
              ccw=[-1 1]*min(AP.ccLagPts,cont2discrete(1*1000/AP.theta(1),osi*.001,'intv',1));
              cci=AP.ccLagPts+[ccw(1):ccw(2)]+1;
              ccp=discrete2cont(cci-AP.ccLagPts,osi*.001,'intv',0);
              hold on
              if length(diag(r(i).thLoeCCMn))>1
                srfc=repmat(nan,length(cci),nAllLFPCh);
                % the lines below take into account non-analyzed channels (which have
                % single nans in the corresponding cells)
                for g=1:nAllLFPCh
                  ix1=min(g,LFPpcInd2);
                  ix2=max(g,LFPpcInd2);
                  if length(r(i).thLoeCCMn{ix1,ix2})>1
                    srfc(1:length(cci),g)=r(i).thLoeCCMn{ix1,ix2}(cci);
                    % symmetry
                    if g>LFPpcInd2, srfc(:,g)=flipud(srfc(:,g)); end
                  else
                    srfc(1:length(cci),g)=nan;
                  end
                end
                % --- contour plot
                [c,cph]=contourf(srfc,nLev);
                clear c;
                axis tight
                set(cph(:),'linestyle','none');
                cph=gca;
                % use the contour plot's x axis as electrode # indicator
                set(cph,'xaxisloc','top','yaxisloc','right','ytick',[],...
                  'xlim',[1 nAllLFPCh],'xtick',[1:nAllLFPCh],'xticklabel',WP.xtl2);
                % x axis limits must be fixed for aligned overlay with surface plot
                ax2=axes('position',get(cph,'position'),'color','none');
                hold on;
                set(gca,'ylim',[ccp(1) ccp(end)],'xlim',[WP.elx(1) WP.elx(end)],...
                  'xtick',WP.elx,'xticklabel',WP.xtl1);
                axis manual;
                tmpMn=cat(1,r(i).thLoeCCPeakTMn{1:LFPpcInd2,LFPpcInd2});
                % symmetry - don't forget to invert
                tmpMn=[tmpMn; -1*cat(1, r(i).thLoeCCPeakTMn{LFPpcInd2,LFPpcInd2+1:end})];
                tmpStd=cat(1, r(i).thLoeCCPeakTStd{1:LFPpcInd2,LFPpcInd2});
                tmpStd=[tmpStd; cat(1, r(i).thLoeCCPeakTStd{LFPpcInd2,LFPpcInd2+1:end})];
                ph=errorbar(WP.elx,tmpMn,tmpStd,[psym '-']);
                set(ph,'color',pcol);
                grid on;
                % line to indicate principal electrode
                lh=line(WP.elx(LFPpcInd2)*[1 1]',get(gca,'ylim'));
                set(lh,'color',pcCol,'linestyle','--','linewidth',4);
                % push to background
                set(gca,'children',circshift(get(gca,'children'),-1))
              end
            end

          case 'gaeCCMn'
            if isfield(r,'gaeCCMn')
              colormap(bone);
              % gamma env (range to plot: the lesser of [+/- 150 ms, AP.ccLagPts])
              if wavefm==1,  
                ccw=[-1 1]*min(AP.ccLagPts,cont2discrete(300,osi*.001,'intv',1));
              else
                ccw=[-1 1]*min(AP.ccLagPts,cont2discrete(150,osi*.001,'intv',1));
              end
              cci=AP.ccLagPts+[ccw(1):ccw(2)]+1;
              ccp=discrete2cont(cci-AP.ccLagPts,osi*.001,'intv',0);
              hold on
              if length(diag(r(i).gaeCCMn))>1
                srfc=repmat(nan,length(cci),nAllLFPCh);
                % the lines below take into account non-analyzed channels (which have
                % single nans in the corresponding cells)
                for g=1:nAllLFPCh
                  ix1=min(g,LFPpcInd2);
                  ix2=max(g,LFPpcInd2);
                  if length(r(i).gaeCCMn{ix1,ix2})>1
                    srfc(1:length(cci),g)=r(i).gaeCCMn{ix1,ix2}(cci);
                    % symmetry
                    if g>LFPpcInd2, srfc(:,g)=flipud(srfc(:,g)); end
                  else
                    srfc(1:length(cci),g)=nan;
                  end
                end
                switch wavefm
                  case 0
                    % --- contour plot
                    [c,cph]=contourf(srfc,nLev);
                    clear c;
                    axis tight
                    set(cph(:),'linestyle','none');
                    cph=gca;
                    % use the contour plot's x axis as electrode # indicator
                    set(cph,'xaxisloc','top','yaxisloc','right','ytick',[],...
                      'xlim',[1 nAllLFPCh],'xtick',[1:nAllLFPCh],'xticklabel',WP.xtl2);
                    % x axis limits must be fixed for aligned overlay with surface plot
                    ax2=axes('position',get(cph,'position'),'color','none');
                    hold on;
                    set(gca,'ylim',[ccp(1) ccp(end)],'xlim',[WP.elx(1) WP.elx(end)],...
                      'xtick',WP.elx,'xticklabel',WP.xtl1);
                    axis manual;
                    tmpMn=cat(1,r(i).gaeCCPeakTMn{1:LFPpcInd2,LFPpcInd2});
                    % symmetry - don't forget to invert
                    tmpMn=[tmpMn; -1*cat(1, r(i).gaeCCPeakTMn{LFPpcInd2,LFPpcInd2+1:end})];
                    tmpStd=cat(1, r(i).gaeCCPeakTStd{1:LFPpcInd2,LFPpcInd2});
                    tmpStd=[tmpStd; cat(1, r(i).gaeCCPeakTStd{LFPpcInd2,LFPpcInd2+1:end})];
                    ph=errorbar(WP.elx,tmpMn,tmpStd,[psym '-']);
                    set(ph,'color',pcol);
                    grid on;
                    % line to indicate principal electrode
                    lh=line(WP.elx(LFPpcInd2)*[1 1]',get(gca,'ylim'));
                    set(lh,'color',pcCol,'linestyle','--','linewidth',4);
                    % push to background
                    set(gca,'children',circshift(get(gca,'children'),-1))
                  otherwise
                    srfc(~isfinite(srfc))=0;
                    [ylim,dy]=pllplot(srfc,'si',osi,'spacing','maxmin');
                    hold on
                    [ylim,dy,yscaleFac,ph]=pllplot(abs(hilbert(srfc)),'spacing','fixed','ylim',ylim,'dy',dy);
                    set(ph,'color','g')
                end
              end
            end

          case 'gaCCMn'
            if isfield(r,'gaCCMn')
              colormap(bone);
              % gamma (range to plot: the lesser of [+/- 1 slow gamma periods, AP.ccLagPts])
              ccw=[-1 1]*min(AP.ccLagPts,cont2discrete(1*1000/AP.gamma(1),osi*.001,'intv',1));
              cci=AP.ccLagPts+[ccw(1):ccw(2)]+1;
              ccp=discrete2cont(cci-AP.ccLagPts,osi*.001,'intv',0);
              hold on
              if length(diag(r(i).gaCCMn))>1
                srfc=repmat(nan,length(cci),nAllLFPCh);
                % the lines below take into account non-analyzed channels (which have
                % single nans in the corresponding cells)
                for g=1:nAllLFPCh
                  ix1=min(g,LFPpcInd2);
                  ix2=max(g,LFPpcInd2);
                  if length(r(i).gaCCMn{ix1,ix2})>1
                    srfc(1:length(cci),g)=r(i).gaCCMn{ix1,ix2}(cci);
                    % symmetry
                    if g>pcIdx, srfc(:,g)=flipud(srfc(:,g)); end
                  else
                    srfc(1:length(cci),g)=nan;
                  end
                end
                switch wavefm
                  case 0
                    % --- contour plot
                    [c,cph]=contourf(srfc,nLev);
                    clear c;
                    axis tight
                    set(cph(:),'linestyle','none');
                    cph=gca;
                    % use the contour plot's x axis as electrode # indicator
                    set(cph,'xaxisloc','top','yaxisloc','right','ytick',[],...
                      'xlim',[1 nAllLFPCh],'xtick',[1:nAllLFPCh],'xticklabel',WP.xtl2);
                    % x axis limits must be fixed for aligned overlay with surface plot
                    ax2=axes('position',get(cph,'position'),'color','none');
                    hold on;
                    set(gca,'ylim',[ccp(1) ccp(end)],'xlim',[WP.elx(1) WP.elx(end)],...
                      'xtick',WP.elx,'xticklabel',WP.xtl1);
                    axis manual;
                    tmpMn=cat(1,r(i).gaCCPeakTMn{1:LFPpcInd2,LFPpcInd2});
                    % symmetry - don't forget to invert
                    tmpMn=[tmpMn; -1*cat(1, r(i).gaCCPeakTMn{LFPpcInd2,LFPpcInd2+1:end})];
                    tmpStd=cat(1, r(i).gaCCPeakTStd{1:LFPpcInd2,LFPpcInd2});
                    tmpStd=[tmpStd; cat(1, r(i).gaCCPeakTStd{LFPpcInd2,LFPpcInd2+1:end})];
                    ph=errorbar(WP.elx,tmpMn,tmpStd,[psym '-']);
                    set(ph,'color',pcol);
                    grid on;
                    % line to indicate principal electrode
                    lh=line(WP.elx(LFPpcInd2)*[1 1]',get(gca,'ylim'));
                    set(lh,'color',pcCol,'linestyle','--','linewidth',4);
                    % push to background
                    set(gca,'children',circshift(get(gca,'children'),-1))
                  otherwise
                    srfc(~isfinite(srfc))=0;
                    pllplot(srfc);
                end

              end
            end

          case 'thCCMn'
            if isfield(r,'thCCMn')
              colormap(bone);
              % III. theta (range to plot: the lesser of [+/- 150 ms, AP.ccLagPts])
              if wavefm==1,  
                ccw=[-1 1]*min(AP.ccLagPts,cont2discrete(300,osi*.001,'intv',1));
              else
                ccw=[-1 1]*min(AP.ccLagPts,cont2discrete(150,osi*.001,'intv',1));
              end
              cci=AP.ccLagPts+[ccw(1):ccw(2)]+1;
              ccp=discrete2cont(cci-AP.ccLagPts,osi*.001,'intv',0);
              hold on
              if length(diag(r(i).thCCMn))>1
                srfc=repmat(nan,length(cci),nAllLFPCh);
                % the lines below take into account non-analyzed channels (which have
                % single nans in the corresponding cells)
                for g=1:nAllLFPCh
                  ix1=min(g,LFPpcInd2);
                  ix2=max(g,LFPpcInd2);
                  if length(r(i).thCCMn{ix1,ix2})>1
                    srfc(1:length(cci),g)=r(i).thCCMn{ix1,ix2}(cci);
                    % symmetry
                    if g>LFPpcInd2, srfc(:,g)=flipud(srfc(:,g)); end
                  else
                    srfc(1:length(cci),g)=nan;
                  end
                end
                switch wavefm
                  case 0
                    % --- contour plot
                    [c,cph]=contourf(srfc,nLev);
                    clear c;
                    axis tight
                    set(cph(:),'linestyle','none');
                    cph=gca;
                    % use the contour plot's x axis as electrode # indicator
                    set(cph,'xaxisloc','top','yaxisloc','right','ytick',[],...
                      'xlim',[1 nAllLFPCh],'xtick',[1:nAllLFPCh],'xticklabel',WP.xtl2);
                    % x axis limits must be fixed for aligned overlay with surface plot
                    ax2=axes('position',get(cph,'position'),'color','none');
                    hold on;
                    set(gca,'ylim',[ccp(1) ccp(end)],'xlim',[WP.elx(1) WP.elx(end)],...
                      'xtick',WP.elx,'xticklabel',WP.xtl1);
                    axis manual;
                    tmpMn=cat(1,r(i).thCCPeakTMn{1:LFPpcInd2,LFPpcInd2});
                    % symmetry - don't forget to invert
                    tmpMn=[tmpMn; -1*cat(1, r(i).thCCPeakTMn{LFPpcInd2,LFPpcInd2+1:end})];
                    tmpStd=cat(1, r(i).thCCPeakTStd{1:LFPpcInd2,LFPpcInd2});
                    tmpStd=[tmpStd; cat(1, r(i).thCCPeakTStd{LFPpcInd2,LFPpcInd2+1:end})];
                    ph=errorbar(WP.elx,tmpMn,tmpStd,[psym '-']);
                    set(ph,'color',pcol);
                    grid on;
                    % line to indicate principal electrode
                    lh=line(WP.elx(LFPpcInd2)*[1 1]',get(gca,'ylim'));
                    set(lh,'color',pcCol,'linestyle','--','linewidth',4);
                    % push to background
                    set(gca,'children',circshift(get(gca,'children'),-1))
                  otherwise
                    srfc(~isfinite(srfc))=0;
                    [ylim,dy]=pllplot(srfc,'si',osi,'spacing','maxmin');
                    hold on
                    [ylim,dy,yscaleFac,ph]=pllplot(abs(hilbert(srfc)),'spacing','fixed','ylim',ylim,'dy',dy);
                    set(ph,'color','g')
                end
              end
            end


          case 'detheCCMn'
            if isfield(r,'detheCCMn')
              colormap(bone);
              % delta-thetaEnv (range to plot: the lesser of [+/- 500 ms, AP.ccLagPts])
              ccw=[-1 1]*min(AP.ccLagPts,cont2discrete(500,osi*.001,'intv',1));
              cci=AP.ccLagPts+[ccw(1):ccw(2)]+1;
              ccp=discrete2cont(cci-AP.ccLagPts,osi*.001,'intv',0);
              % limits [-.33 .27]
              v=((0:nLev-1)/(nLev-1)-.55)*.1*6;
              hold on
              if ~isempty(r(i).detheCCMn)
                srfc=repmat(nan,length(cci),nAllLFPCh);
                srfc(:,LFPccInd)=r(i).detheCCMn(cci,:);
                % --- contour plot
                [c,cph]=contourf(srfc,v);
                clear c;
                axis tight
                set(cph(:),'linestyle','none');
                set(gca,'clim',v([1 end]));
                cph=gca;
                % use the contour plot's x axis as electrode # indicator
                set(cph,'xaxisloc','top','yaxisloc','right','ytick',[],...
                  'xlim',[1 nAllLFPCh],'xtick',[1:nAllLFPCh],'xticklabel',WP.xtl2);
                % x axis limits must be fixed for aligned overlay with surface plot
                ax2=axes('position',get(cph,'position'),'color','none');
                hold on;
                set(gca,'ylim',[ccp(1) ccp(end)],'xlim',[WP.elx(1) WP.elx(end)],...
                  'xtick',WP.elx,'xticklabel',WP.xtl1);
                axis manual;
                tmpMn=tempo;
                tmpStd=tempo;
                tmpMn(LFPccInd)=r(i).detheCCPeakTMn;
                tmpStd(LFPccInd)=r(i).detheCCPeakTStd;
                ph=errorbar(WP.elx,tmpMn,tmpStd,[psym '-']);
                set(ph,'color',pcol);
                grid on;
                % line to indicate principal electrode
                lh=line(WP.elx(LFPpcInd2)*[1 1]',get(gca,'ylim'));
                set(lh,'color',pcCol,'linestyle','--','linewidth',4);
                % push to background
                set(gca,'children',circshift(get(gca,'children'),-1))
              end
            end



          case 'deCCMn'
            if isfield(r,'deCCMn')
              colormap(bone);
              % delta (range to plot: the lesser of [+/- 500 ms, AP.ccLagPts])
              ccw=[-1 1]*min(AP.ccLagPts,cont2discrete(500,osi*.001,'intv',1));
              cci=AP.ccLagPts+[ccw(1):ccw(2)]+1;
              ccp=discrete2cont(cci-AP.ccLagPts,osi*.001,'intv',0);
              hold on
              if length(diag(r(i).deCCMn))>1
                srfc=repmat(nan,length(cci),nAllLFPCh);
                % the lines below take into account non-analyzed channels (which have
                % single nans in the corresponding cells)
                for g=1:nAllLFPCh
                  ix1=min(g,LFPpcInd2);
                  ix2=max(g,LFPpcInd2);
                  if length(r(i).deCCMn{ix1,ix2})>1
                    srfc(1:length(cci),g)=r(i).deCCMn{ix1,ix2}(cci);
                    % symmetry
                    if g>LFPpcInd2, srfc(:,g)=flipud(srfc(:,g)); end
                  else
                    srfc(1:length(cci),g)=nan;
                  end
                end
                % --- contour plot
                [c,cph]=contourf(srfc,nLev);
                clear c;
                axis tight
                set(cph(:),'linestyle','none');
                cph=gca;
                % use the contour plot's x axis as electrode # indicator
                set(cph,'xaxisloc','top','yaxisloc','right','ytick',[],...
                  'xlim',[1 nAllLFPCh],'xtick',[1:nAllLFPCh],'xticklabel',WP.xtl2);
                % x axis limits must be fixed for aligned overlay with surface plot
                ax2=axes('position',get(cph,'position'),'color','none');
                hold on;
                set(gca,'ylim',[ccp(1) ccp(end)],'xlim',[WP.elx(1) WP.elx(end)],...
                  'xtick',WP.elx,'xticklabel',WP.xtl1);
                axis manual;
                tmpMn=cat(1,r(i).deCCPeakTMn{1:LFPpcInd2,LFPpcInd2});
                % symmetry - don't forget to invert
                tmpMn=[tmpMn; -1*cat(1, r(i).deCCPeakTMn{LFPpcInd2,LFPpcInd2+1:end})];
                tmpStd=cat(1, r(i).deCCPeakTStd{1:LFPpcInd2,LFPpcInd2});
                tmpStd=[tmpStd; cat(1, r(i).deCCPeakTStd{LFPpcInd2,LFPpcInd2+1:end})];
                ph=errorbar(WP.elx,tmpMn,tmpStd,[psym '-']);
                set(ph,'color',pcol);
                grid on;
                % line to indicate principal electrode
                lh=line(WP.elx(LFPpcInd2)*[1 1]',get(gca,'ylim'));
                set(lh,'color',pcCol,'linestyle','--','linewidth',4);
                % push to background
                set(gca,'children',circshift(get(gca,'children'),-1))
              end
            end

          case 'thgaeCCMn'
            if isfield(r,'thgaeCCMn')
              colormap(bone);
              % theta-gammaEnv (range to plot: [+/- 150 ms])
              ccw=[-1 1]*min(AP.ccLagPts,cont2discrete(150,osi*.001,'intv',1));
              cci=AP.ccLagPts+[ccw(1):ccw(2)]+1;
              ccp=discrete2cont(cci-AP.ccLagPts,osi*.001,'intv',0);
              hold on
              if ~isempty(r(i).thgaeCCMn)
                srfc=repmat(nan,length(cci),nAllLFPCh);
                srfc(:,LFPccInd)=r(i).thgaeCCMn(cci,:);
                switch wavefm
                  case 0
                    % --- contour plot
                    [c,cph]=contourf(srfc,nLev);
                    clear c;
                    axis tight
                    set(cph(:),'linestyle','none');
                    cph=gca;
                    % use the contour plot's x axis as electrode # indicator
                    set(cph,'xaxisloc','top','yaxisloc','right','ytick',[],...
                      'xlim',[1 nAllLFPCh],'xtick',[1:nAllLFPCh],'xticklabel',WP.xtl2);
                    % x axis limits must be fixed for aligned overlay with surface plot
                    ax2=axes('position',get(cph,'position'),'color','none');
                    hold on;
                    set(gca,'ylim',[ccp(1) ccp(end)],'xlim',[WP.elx(1) WP.elx(end)],...
                      'xtick',WP.elx,'xticklabel',WP.xtl1);
                    axis manual;
                    tmpMn=tempo;
                    tmpStd=tempo;
                    tmpMn(LFPccInd)=r(i).thgaeCCPeakTMn;
                    tmpStd(LFPccInd)=r(i).thgaeCCPeakTStd;
                    ph=errorbar(WP.elx,tmpMn,tmpStd,[psym '-']);
                    set(ph,'color',pcol);
                    grid on;
                    % line to indicate principal electrode
                    lh=line(WP.elx(LFPpcInd2)*[1 1]',get(gca,'ylim'));
                    set(lh,'color',pcCol,'linestyle','--','linewidth',4);
                    % push to background
                    set(gca,'children',circshift(get(gca,'children'),-1))
                  otherwise
                    srfc(~isfinite(srfc))=0;
                    pllplot(srfc);
                end
              end
            end % if:isfield
          case {'rawDePEMn','rawThPEMn','rawThNarrowPEMn','rawBePEMn','rawGaPEMn','rawRiPEMn'}
            ylab='P (mV^2)';
            % assume that if theta power exists power had also been calculated in other bands
            if length(diag(r(i).rawThNarrowPEMn))>1
              eval(['tmpMn=cat(2,r(i).' rv{rvi} '{AP.dixie});'])
              eval(['tmpStd=cat(2,r(i).' rv{rvi}(1:end-2) 'Std{AP.dixie});'])
              ph=errorbar(WP.elx,tmpMn,zeros(size(tmpMn)),tmpStd,[psym '-']);
              set(ph,'color',pcol);
            end
            % line (below) will not be drawn if lower y limit is <0 and then yscale set to
            % log, so set y axis limits manually
            nicexyax;
            % yl=get(gca,'ylim');
            % yl=yl.*[.8 1.2];
            % set(gca,'ylim',yl,'yscale','log','xtick',WP.elx,'xticklabel',WP.xtl1);
            set(gca,'xtick',WP.elx,'xticklabel',WP.xtl1);
            ylabel(ylab);
            % line to indicate principal electrode
            lh=line(WP.elx(LFPpcInd2)*[1 1]',get(gca,'ylim'));
            set(lh,'color',pcCol,'linestyle',':');
            % push to background
            set(gca,'children',circshift(get(gca,'children'),-1))
            % second axis for electrode #
            ax2=axes('position',get(gca,'position'),'color','none','xaxisloc','top',...
              'yaxisloc','right','ylim',get(gca,'ylim'),'ytick',get(gca,'ytick'),...
              'xlim',get(gca,'xlim'),'xtick',WP.elx,'xticklabel',WP.xtl2,'yticklabel',[]);
        end % switch:rv

        th=ultext([DS.aName '; ' DS.abfFn],.05,'color','b','fontweight','b','fontsize',8);
        % once more, tiny offset
        th=ultext([DS.aName '; ' DS.abfFn],.053,'color',[.7 .7 .1],'fontweight','b','fontsize',8);
        drawnow
        if ci==n2 & ri==n1
          % Pump 'em up, Arnie
          tmp1=get(gcf,'children');
          for pipi=1:length(tmp1)
            if ~strcmpi(get(tmp1(pipi),'tag'),'colorbar')
              rexy('ax',tmp1(pipi),'xfac',1.1,'yfac',1.15);
            end
          end
          % print?
          if ~isempty(printas{1}),
            for pipi=1:length(printas)
              pa=printas{pipi};
              if strfind(pa,'ps'), ext='.ps';
              elseif strfind(pa,'jpeg'), ext='.jpg';
              else ext='';
              end
              try
                print(pa,[WP.rootPath figdir '\' figName ext]);
              catch
                print(pa,[WP.rootPath '\' figName ext]);
              end
            end
          end

        end % if:last file reached
      end % for:nRv
    end % for:behav
  end % for:rows of ANPAR=concs
end % for:cols of ANPAR=experiments


% ------------ local func ----------------- local func --------------------
% ------------ local func ----------------- local func --------------------

function figha=mkfig(ftag)
figha=findobj('tag',ftag);
if isempty(figha), figha=figure;
else  figure(figha);
end
tmpScrSz=get(0,'Screensize')*.7;
tmpScrSz([1 2])=tmpScrSz([1 2])+tmpScrSz([3 4])*.10+25*rand;
set(figha,'position',tmpScrSz,'tag',ftag,'name',ftag,...
  'color',[.9 .9 1],'numbertitle','off');
