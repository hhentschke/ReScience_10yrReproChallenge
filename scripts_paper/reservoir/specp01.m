function specp01
% generates plots of spectra (ctrl and drug in one plot; two channels and
% two behaviors)

global DS AP

% load data from original results file (excruciatingly slow) and save extracted data
saveMd='original';  
% load previously extracted data (fast)
% saveMd='extracted';  

rv={'rawCohMn'};
rv={'rawPMn'};
rvi=1;
      
behav={'immobile','exploring'};
locChan={'IN 6';'IN 11'};
nLocChan=length(locChan);
% x axis limits
xl=[1 100];

pcol={[.6 .6 .6],'k'};

% --- prepare graphics
labelscale('fontSz',8,'scaleFac',.5,'lineW',1.5,'markSz',8); 
ornt='portrait';
figdir='c:\projects\rmouse\paper_atropine\rawFig\';
printas='-dpsc2';[];
figName=mfilename;

figure(1), clf, orient(ornt); 
% ------ main loop: files
for fi=1:2
  % --- paths, channels
  rmouse_ini;
  cd([WP.rootPath '\beta3_wtko\wt0001_04708']);
  switch fi
    case 1
      a001_r1;
    case 2
      % a003_r1;    
      a002_exc1;    
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
          case {'rawPMn'}
            tmpc=cat(2,tmpc{AP.dixie(inter_ix)});      
            % cut down
            tmpc=tmpc(cix,:);
            c=cat(3,c,tmpc);
            ylab='PSD (mV^2/Hz)';
          case {'rawCohMn'}      
            error('not yet done');
            % the trick: redundant combinations are empty
            tmpc=cat(1,c{setdiff(inter_ix,AP.LFPpcInd2),AP.LFPpcInd2});
            tmpc=cat(1,tmpc,c{AP.LFPpcInd2,inter_ix});
            tmpc=reshape(tmpc,[length(F),nLocChan]);
            % cut down
            tmpc=tmpc(cix,:);
            c=tmpc;
            ylab='Coherence';
        end
      end
      save([figdir DS.abfFn '_tmpDat'],'ix1','ix2','xax','cix','F','c','ylab');
      
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
      subplot(nLocChan,length(behav),(loci-1)*length(behav)+bi); 
      hold on
    
      ph=plot(xax,c(:,loci,bi));
      set(ph,'color',pcol{fi});
      if strfind(rv{rvi},'Coh')
        set(gca,'xscale','log','ylim',[0 1.1]);
      else
        set(gca,'xscale','log','yscale','log');
        % set(gca,'xscale','log');        
      end
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

% for each channel, set y axes to same value
chld=get(gcf,'children');

mima=[get(chld(1),'ylim')  get(chld(2),'ylim')];
mima=[min(mima) max(mima)];
set(chld(1),'ylim',mima);
set(chld(2),'ylim',mima);

mima=[get(chld(3),'ylim')  get(chld(4),'ylim')];
mima=[min(mima) max(mima)];
set(chld(3),'ylim',mima);
set(chld(4),'ylim',mima);

% only now can we set the rectangles defining frequency bands
for ci=1:length(chld)
  subplot(chld(ci))
  
  ylim=get(gca,'ylim');
  ph=patch([AP.delta+[0 -.1] fliplr(AP.delta+[0 -.1] )],reshape([ylim;ylim],1,4),[.8 .8 .8]);
  set(ph,'linestyle','none')
  % push to background
  set(gca,'children',circshift(get(gca,'children'),-1))
  
  ylim=get(gca,'ylim');
  ph=patch([AP.theta fliplr(AP.theta)],reshape([ylim;ylim],1,4),[.8 .8 .8]);
  set(ph,'linestyle','none')
  % push to background
  set(gca,'children',circshift(get(gca,'children'),-1))
  
  ph=patch([AP.beta fliplr(AP.beta)],reshape([ylim;ylim],1,4),[.8 .8 .8]);
  set(ph,'linestyle','none')
  % push to background
  set(gca,'children',circshift(get(gca,'children'),-1))
  
  ph=patch([AP.gamma fliplr(AP.gamma)],reshape([ylim;ylim],1,4),[.8 .8 .8]);
  set(ph,'linestyle','none')
  % push to background
  set(gca,'children',circshift(get(gca,'children'),-1))
end

if ~isempty(printas), 
  print(printas,[figdir figName '_' rv{rvi}]); 
end




