function ccmovie(r,varargin)
% ** function ccmovie(r,varargin)
% 'CC movie' generating function working on 3D variables as generated 
% by function catf.m
%                         >>> INPUT VARIABLES >>>
% NAME             TYPE/DEFAULT          DESCRIPTION
% r                struct array          results struct of catf.m
% rv               cell array of chars   any (combination) of fields of r containing
%                                        segment-wise computed parameters, like 'thCCPeak'
% chanComb         char arr, 'neigh'     'neigh' - nearest neighbor (line plots)
%                                        'princ' - all vs. principal channel (line plots)
%                                        'all'   - all, color-coded matrix plot
%                     
%                         <<< OUTPUT VARIABLES <<<
% NAME             TYPE/DEFAULT           DESCRIPTION



% extension: look for recurring patterns using PC, make movie
% more than one variable at once
% options: interpolation, smoothing
% show time, behavior

% clear separation of tasks:
% - pattern recognition routine: e.g. PC computing routine
% - option in cc movie: plot (i) one state variable as running time plot 
% (ii) two state variables as running scatter plot


global DS AP 

chanComb='all';
pvpmod(varargin);

% --- preliminaries: local copies of AP and DS & checks
AP=r.AP;
DS=r.DS;
rmouse_APcheck;
rawCh=rmouse_chan;
nCh=nLFPCh;


% --- preliminaries: initialize figure
ftag='ccmov';
fh=findobj('tag',ftag);
if isempty(fh), fh=figure;
else  figure(fh);
end
tmpScrSz=get(0,'Screensize');
tmpScrSz=tmpScrSz*.8;
tmpScrSz(1)=tmpScrSz(1)+tmpScrSz(3)*.1;
tmpScrSz(2)=tmpScrSz(2)+tmpScrSz(4)*.1;  
set(fh,'position',tmpScrSz,...
  'tag',ftag,...
  'name',ftag,...
  'color',[0.27 0.27 .4],...
  'numbertitle','off');
clf;
orient landscape;

% --- preliminaries: graphics settings
labelscale('fontSz',7,'scaleFac',1.0,'lineW',1.0,'markSz',4); 
% set up compound colormap
nColors=128;
colormap([jet(nColors); coma('redblue','ncols',nColors)]);
% # of rows = # of CC types
nRows=length(rv);
% # of cols = 2 (one for amplitude, one for lag)
nCols=2;
xmarg=.04;
ymarg=.035;
xpos=(0:nCols-1)*((1-2*xmarg)/nCols)+2*xmarg;
xlen=(1-2*xmarg)/nCols-2*xmarg;
% leave additional space (.05)for title subplot
ypos=fliplr((0:nRows-1)*((1-.05-2*ymarg)/nRows)+2*ymarg);
ylen=(1-.05-2*ymarg)/nRows-2*ymarg;
ccprop={'ampl','lag'};
cblabFstr={'%1.1f','%3.0f'};
% title 
% in title string replace underscores (tex interpreter (?) makes text to subscript)
tmpfn1=DS.abfFn; tmpfn1(strfind(tmpfn1,'_'))='-'; 
tmpfn2=AP.bScoreFn; tmpfn2(strfind(tmpfn2,'_'))='-';         
figTitle=[DS.aName ' (' tmpfn1 ' + ' tmpfn2 ')'];
subplot('position',[xmarg .95 1-2*xmarg .02]);
set(gca,'xlim',[0 1],'ylim',[0 1]);
axis off
th=text(.5,.5,figTitle,'fontsize',12,'fontweight','bold','horizontalalignment','center');

% --- preliminaries: bounds for parameters
for sti=1:nRows
  switch rv{sti}
    case {'deCCPeak','deCCPeakT'}
      maxLag=250;
      minAmp=.0;
    case {'thCCPeak','thCCPeakT'}
      maxLag=60;
      minAmp=.0;
    case {'thLoeCCPeak','thLoeCCPeakT','thHieCCPeak','thHieCCPeakT'}
      maxLag=180;
      minAmp=.0;
    case {'gaCCPeak','gaCCPeakT'}
      maxLag=10;
      minAmp=.0;
    case {'gaeCCPeak','gaeCCPeakT'}
      maxLag=10;
      minAmp=.0;
    otherwise
      error('illegal streamType in CC matrix plot');
  end
end
% --- load variables & put in array
switch chanComb
  case 'neigh'
  case 'all'
    % the big scaling: since a colormap applies to a whole figure (and not a
    % single axis) we have to use a compound color map and scale the data such that
    % CC lags and amplitudes are converted to sets of nonoverlapping indices
    % 1. amplitudes: map [0 1] to [1 nColors]
    cFac(1)=nColors-1;
    cOffs(1)=1;
    % 2. lags: map [-maxLag maxLag] to [nColors+1 2*nColors]
    cFac(2)=(nColors-1)/(2*maxLag);
    cOffs(2)=cFac(1)+cOffs(1)+nColors/2;
    
    % --- pruning, plucking, curbing
    % 1. amplitudes
    % CC peak values should be within [0 1], but to prevent indexing errors make
    % sure that is really so:
    ccAmp(ccAmp>1)=1;
    ccAmp(ccAmp<0)=0;
    ccAmp=round(cOffs(1)+cFac(1)*ccAmp);
    % 2. lags:
    % same story: deal with outliers (set to nan here)
    ccLag(abs(ccLag)>maxLag)=nan;
    ccLag=round(cOffs(2)+cFac(2)*ccLag);
end


switch chanComb
  case 'neigh'
  case 'all'
    % set up subplots
    ct=0;
    for sti=1:nRows
      for ii=1:2
        ct=ct+1;
        sph(ct)=subplot('position',[xpos(1+mod(ii-1,2)) ypos(sti) xlen ylen]);
        if ii==1
          ih(ct)=image(ccAmp(:,:,1));
        else
          ih(ct)=image(ccLag(:,:,1));
        end
        set(ih,'CDataMapping','direct')
        set(gca,'xaxisloc','top',...
          'ylim',[.5 nAllLFPCh+.5],'ytick',[1:nAllLFPCh],'yticklabel',WP.xtl2,...
          'xlim',[.5 nAllLFPCh+.5],'xtick',[1:nAllLFPCh],'xticklabel',WP.xtl2);
        axis manual
        axis square
        % inflate plots
        %rexy('xfac',1.4,'yfac',1.4);
        title([r(rix).segmentType ';' strmType ';' ccprop{ii}]);
      end
    end
    
    for k=1:r(rix).ni
      ct=0;
      for sti=1:nRows
        for ii=1:2
          ct=ct+1;
          subplot(sph(ct));
          if ii==1
            set(ih(ct),'CData',ccAmp(:,:,k));
          else
            set(ih(ct),'CData',ccLag(:,:,k));
          end
          % lines to indicate principal electrode
          lh=line(AP.LFPpcInd2*[1 1]',[.5 AP.LFPpcInd2]);
          set(lh,'color',[.6 .6 .6],'linestyle','-');
          lh=line([AP.LFPpcInd2 nAllLFPCh+.5],AP.LFPpcInd2*[1 1]');
          set(lh,'color',[.6 .6 .6],'linestyle','-');
          % in last frame put colorbar
          if k==r(rix).ni
            cbh=colorbar('vert');
            set(cbh,'ylim',[(ii-1)*nColors  ii*nColors]);
            % set ticks/labels at/to meaningful positions/values
            set(cbh,'ytick',[ii-1 ii-.5  ii]*nColors,...
              'yticklabel',num2str(((([ii-1 ii-.5  ii]*nColors)-cOffs(ii))/cFac(ii))',cblabFstr{ii}));
          end
        end
      end
      drawnow;
      pause(.3)
    end
end


% getframe
% movie