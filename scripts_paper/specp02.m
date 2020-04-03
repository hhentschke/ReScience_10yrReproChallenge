function specp01
% plots raw and gammaEnv spectra (averages as computed by rmouse) in one plot.
% One behavior, one channel, ctrl

global DS AP

% load data from original results file (excruciatingly slow) and save extracted data
saveMd='original';
% load previously extracted data (fast)
saveMd='extracted';


rv={'rawPMn','gaePMn'};

behav={'exploring'};
locChan={'IN 11'};
nLocChan=length(locChan);
% x axis limits
xl=[2 90];

% --- prepare graphics
labelscale('fontSz',8,'scaleFac',.25,'lineW',1.5,'markSz',8);
ornt='portrait';
figdir='d:\projects\hippocampus\rmouse\paper_atropine\rawFig\';
printas='-dpsc2';
% printas=[];

figName=mfilename;

figure(1), clf, orient(ornt);

fi=1;
% --- paths, channels
rmouse_ini;
cd([WP.rootPath '\beta3_wtko\wt0001_04708']);
a001_r1;
AP_beta3_wtko;
rmouse_apcheck;
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
    gaeF=r(1).gaeF;
    % --- all about time: axes
    xax=F;
    gaexax=gaeF;
    if ~isempty(xl)
      cix=find(xax>=xl(1) & xax<=xl(2));
      xax=xax(cix);
      gaecix=find(gaexax>=xl(1) & gaexax<=xl(2));
      gaexax=gaexax(gaecix);
    else
      cix=1:length(xax);
      gaecix=1:length(gaexax);
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
    for rvi=1:length(rv)
      c=[];
      for bi=sort(rix)
        eval(['tmpc=r(bi).' rv{rvi} ';']);
        % structure of variables: freq|channel|behav
        tmpc=cat(2,tmpc{AP.dixie(inter_ix)});
        switch rv{rvi}
          case {'rawPMn'}
            % cut down
            tmpc=tmpc(cix,:);
            c=cat(3,c,tmpc);
          case {'rawPMn','gaePMn'}
            % cut down
            tmpc=tmpc(gaecix,:);
            c=cat(3,c,tmpc);
        end
      end
      eval([rv{rvi} '=c;']);
      ylab='PSD (mV^2/Hz)';
      if rvi==1
        save([figdir DS.abfFn '_tmpDat'],'ix1','ix2','xax','cix','F','gaexax','gaecix','gaeF',rv{rvi},'ylab');
      else
        save([figdir DS.abfFn '_tmpDat'],'ix1','ix2','xax','cix','F','gaexax','gaecix','gaeF',rv{rvi},'ylab','-append');
      end
    end

  case 'extracted'
    load([figdir DS.abfFn '_tmpDat']);
end

% --- plot
for loci=1:nLocChan
  for h=1:length(ix1)
    % this substitutes 60 Hz peaks and harmonics by mean of neighboring bins
    rawPMn(ix1{h},:)=repmat(mean(rawPMn(ix2{h},:),1),[3 1]);
  end
  for bi=1:length(behav)
    hold on

    ph=plot(xax,rawPMn(:,loci,bi),'k');
    ph=plot(gaexax,gaePMn(:,loci,bi),'k');

    set(gca,'xscale','log','yscale','log');
    % set(gca,'xscale','log');

    set(gca,'xlim',[max(xax(1),xl(1)) min(xax(end),xl(end))]);
    % niceyax etc does not work with log scaling
    ymin=min(gaePMn(:,loci,bi));
    ymax=max(rawPMn(:,loci,bi));
    ylim=[ymin ymax+(ymax-ymin)*.2];
    set(gca,'ylim',ylim)
    set(gca,'xtick',[1 2.5 5 10 20 40 80]);
    set(gca,'ytick',10.^[-6:0]);
    
    
    ylim=get(gca,'ylim');
    ph=patch([AP.theta fliplr(AP.theta)],reshape([ylim;ylim],1,4),[.8 .8 .8]);
    set(ph,'linestyle','none')
    % push to background
    set(gca,'children',circshift(get(gca,'children'),-1))

    ph=patch([AP.gamma fliplr(AP.gamma)],reshape([ylim;ylim],1,4),[.8 .8 .8]);
    set(ph,'linestyle','none')
    % push to background
    set(gca,'children',circshift(get(gca,'children'),-1))


  end % for:behaviors
end % for:channels

if ~isempty(printas),
  print(printas,[figdir figName]);
end




