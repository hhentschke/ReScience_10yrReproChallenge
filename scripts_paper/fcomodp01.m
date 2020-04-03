function fComodp01
% makes a plot of frequency comodulograms (ctrl and drug in one separate figures; 
% principal channel and exploring only)

global DS AP

% load data from original results file (slow) and save extracted data
saveMd='original';  
% load previously extracted data (fast)
saveMd='extracted';  


rv={'fComod'};
rvi=1;
      
behav={'exploring'};
locChan={'IN 11'};
nLocChan=length(locChan);
% axis limits
xl=[0 100];


% --- prepare graphics
labelscale('fontSz',8,'scaleFac',.6,'lineW',.5,'markSz',8); 


ornt='portrait';
figdir='c:\projects\hippocampus\rmouse\paper_atropine\rawFig\';
printas='-djpeg99';[];
printas='-dpsc2';
% printas=[];

% ------ main loop: files
cnt=0;
for fi=1:2
  % --- paths, channels
  rmouse_ini;
  cd([WP.rootPath '\beta3_wtko\wt0001_04708']);
  switch fi
    case 1
      a001_r1;
    case 2
      a003_r1;    
  end
  AP_beta3_wtko;
  rmouse_APcheck;
  rawCh=rmouse_chan;
  % locate desired channels within results structure
  for i=1:nLocChan
    % indices for intra-electrode CC (AP.LFPInd is not needed..)
    intra_ix(i)=strmatch(locChan{i},AP.rawChAnNm(AP.LFPInd),'exact');
    % indices for inter-electrode CC
    inter_ix(i)=strmatch(locChan{i},DS.rawCh(AP.allLFPIdx,1),'exact');  
    % principal channels should be clear
  end
  if isempty(strfind(DS.dpath,':')), DS.dpath=[WP.rootPath DS.dpath]; end
  
  % --- load 
  switch saveMd
    case 'original'
      load([AP.resPath '\' AP.resFn],'r');
      F=r(1).comF;
      % --- all about time: axes
      xax=F;
      if ~isempty(xl)
        cix=find(xax>=xl(1) & xax<=xl(2));
        xax=xax(cix);
      else
        cix=1:length(xax);
      end
      
      % --- extract 
      % determine index to results obtained with specified behavior in r
      % note that intersect will sort alphabetically, which we dont want)
      [nada,rix]=intersect({r(:).segmentType},behav);
      c=[];
      for bi=sort(rix)
        % structure of c: freq|channel|behav
        eval(['tmpc=r(bi).' rv{rvi} ';']);
        switch rv{rvi}
          case {'fComod'}
            c=tmpc(:,:,intra_ix);
        end
      end
      save([figdir DS.abfFn '_tmpDat'],'xax','cix','F','c');
      clear tmpc r
      
    case 'extracted'
      load([figdir DS.abfFn '_tmpDat']);
  end      
  
  % --- plot
  
  % assume that bulk of correlation coeffs remain in that interval
  ccLim=[-.20 .35];
  nCol=255;

  for loci=1:nLocChan
    for bi=1:length(behav)
      cnt=cnt+1;
      figName=[mfilename '_' DS.abfFn '_' behav{bi} '_' locChan{loci}];
      figure(cnt), clf;
      orient(ornt);
      
      if 0
        % colormap(coma('turquoisamber','ncols',nCol));
        % colormap(coma('blueblackorange','ncols',nCol));
        % colormap(coma('mib','ncols',nCol));
        cm=colormap;
      else
        cm=[];
        % colormap
        partCm=colormap('hot');
        % get rid of white parts
        partCm=partCm(1:round(size(partCm,1)*.85),:);
        % interpolate to get number of colors right
        [n1,n2]=size(partCm);
        for g=n2:-1:1
          cm(:,g)=interp1((1:n1)',partCm(:,g),linspace(1,n1,nCol/2)');
        end
        clear partCm
        % swap red and blue column, then flipud & concatenate
        cm=[flipud((cm(:,[3 2 1])).^.35); cm];
      end

      % trim colormap such that 0 corr is black
      if abs(ccLim(1))<ccLim(2)
        cutIx=floor((1-abs(ccLim(1))/ccLim(2))*(nCol/2));
        cm=cm(cutIx:end,:);
        colormap(cm);
      end
      
      pc=tril(c(cix,cix,loci),-1);
      if 0
        % set upper triangular matrix to mean value for symmetric colormaps
        % with white in middle
        pc(~pc)=mean(ccLim);
      else
        % curb max values and set upper triangular matrix to any high value 
        % for symmetric colormaps with black in middle..
        pc(pc>ccLim(2))=ccLim(2);
        pc(~pc)=ccLim(2)+diff(ccLim);
        % ..append white as last color..
        cm=colormap; 
        colormap([cm; 1 1 1]);
        % and set clim such that high values in upper tringle appear white
        ccLim(2)=ccLim(2)+diff(ccLim)/(size(cm,1)-1);
      end
      if 1
        cph=imagesc(F(cix),F(cix),pc,ccLim);
      else
        cph=pcolor(F(cix),F(cix),pc);
        set(cph,'Linestyle','none')
        caxis(ccLim);
        set(gca,'yscale','log','xscale','log','ydir','reverse');
      end
      axis square
      set(gca,'xtick',[10:20:F(end)],'ytick',[10:20:F(end)],...
        'xaxislocation','bott');
      grid on
      colorbar;
      yl=F(cix([1 end]));
      ph=patch([8;8;11;11],[40;yl(2)-.5;yl(2)-.5;40],'r');
      set(ph,'linewidth',1','facecolor','none','edgecolor',[.8 .8 .8],'linestyle',':');
      
      rexy('ax','gca','xfac',.5,'yfac',.5);

      if ~isempty(printas),
        print(printas,[figdir figName '_' rv{rvi}]);
      end

    end % for:behaviors
  end % for:channels
end % for:files






