function specp03
% generates plots of spectra (ctrl, iso and iso+atropine in one plot; one channel and
% two behaviors for control case)

global DS AP

% load data from original results file (excruciatingly slow) and save extracted data
saveMd='original';  
% load previously extracted data (fast)
% saveMd='extracted';  

rv={'gaePMn'};
rv={'rawPMn'};
rvi=1;
      
behav={'immobile','exploring'};
% behav={'immobile'};

locChan={'IN 3'};
nLocChan=length(locChan);
% x axis limits
xl=[2 80.1];

pcol={[.6 .6 .6],'k','compound'};

% --- prepare graphics
labelscale('fontSz',8,'scaleFac',.75,'lineW',1.5,'markSz',8); 
ornt='portrait';
figdir='d:\projects\rmouse\paper_atropine\rawFig\';
printas='-dpsc2';[];
figName=[mfilename '_' rv{rvi}];

figure(1), clf, orient(ornt); 
% array of plot handles
pha=[];

% ------ main loop: files
for fi=1:3
  % --- paths, channels
  rmouse_ini;
  cd([WP.rootPath '\beta3_wtko\isoatr\wt2141_02523']);
  switch fi
    case 1
      a000_r1;
    case 2
      a001_r1;
    case 3
      a002_r1;    
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
      F=r(1).F;
      % --- all about time: axes
      xax=F;
      if ~isempty(xl)
        cix=find(xax>=xl(1) & xax<=xl(2));
        xax=xax(cix);
      else
        cix=1:length(xax);
      end
      
      % preparations for cutting out line hum
      humF=60:60:xax(end);
      for h=1:length(humF)
        % indices to immediately adjacent bins
        [nix,ix]=min(abs(xax-humF(h)));
        % the ones to be replaced
        ix1a=ix-1:ix+1;
        ix1a(ix1a<1 | ix1a>length(xax))=[];
        ix1{h}=ix1a;
        % the ones to compute mean from & replace with 
        ix2a=[ix-5:ix-2 ix+2:ix+5];
        ix2a(ix2a<1 | ix2a>length(xax))=[];            
        ix2{h}=ix2a;
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
          case {'rawPMn','gaePMn'}
            if (iscell(tmpc) && isempty(tmpc{1})) || isempty(tmpc)
              tmpc=repmat(nan,length(cix),1);
            else
              tmpc=cat(2,tmpc{AP.dixie(inter_ix)});
              % cut down
              tmpc=tmpc(cix,:);
            end
            c=cat(3,c,tmpc);
            ylab='PSD (mV^2/Hz)';
        end
      end
      save([figdir DS.abfFn rv{rvi} '_tmpDat'],'ix1','ix2','xax','cix','F','c','ylab');
    case 'extracted'
      load([figdir DS.abfFn '_tmpDat']);
  end      
  
  % --- plot
  for loci=1:nLocChan
    for h=1:length(ix1)
      % this substitutes 60 Hz peaks and harmonics by mean of neighboring bins
      c(ix1{h},:)=repmat(mean(c(ix2{h},:),1),[3 1]);          
    end
    for bi=1:length(behav)
      subplot(2,length(behav),(loci-1)*length(behav)+bi); 
      hold on
      ph=plot(xax,c(:,loci,bi));
      if bi==1, pha(fi)=ph; end
      if strcmpi(pcol{fi},'compound')
        set(ph,'color','k');
        % overlay thin trace in white
        ph=plot(xax,c(:,loci,bi));
        set(ph,'linewidth',get(ph,'linewidth')*.6,'color','w');
      else
        set(ph,'color',pcol{fi});
      end
      set(gca,'xscale','log','yscale','log');
      set(gca,'xlim',xax([1 end]));    
      % niceyax etc does not work with log scaling
      ymin=min(c(:,loci,bi));
      ymax=max(c(:,loci,bi));
      if fi==1
        set(gca,'ylim',[ymin ymax+(ymax-ymin)*.2])
      else
        ylim=get(gca,'ylim');
        if ymin<ylim(1)
          set(gca,'ylim',[ymin ylim(2)]);
        end
        if ymax>ylim(2)
          set(gca,'ylim',[ylim(1) ymax+(ymax-ymin)*.2]);
        end
%         ylim=get(gca,'ylim');
%         set(gca,'ylim',[ylim(1) ylim(2)+diff(ylim)*.2]);
        set(gca,'xtick',[1 2.5 5 10 20 40 80]);
        set(gca,'ytick',10.^[-6:0])
%         xlabel('Freq (Hz)');
%         ylabel(ylab);
      end
      % rexy('ax',gca,'xfac',.5,'yfac',.5);
    end % for:behaviors
  end % for:channels
end % for:files

% for each subplot set y axes to same value
chld=get(gcf,'children');

mima=[get(chld(1),'ylim')  get(chld(2),'ylim')];
mima=[min(mima) max(mima)];
set(chld(1),'ylim',mima);
set(chld(2),'ylim',mima);

sph=subplot(2,2,3)
legend(sph,pha,{'ctrl','iso','iso+atr'});

if ~isempty(printas), 
  print(printas,[figdir figName '_' rv{rvi}]); 
end




