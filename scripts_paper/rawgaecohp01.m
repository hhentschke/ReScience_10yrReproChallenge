function rawgaecohp01
% plots coherence of raw and gammaEnv in one plot.
% One behavior, ctrl

global DS AP

% load data from original results file (slow) and save extracted data
saveMd='original';
% load previously extracted data (fast)
saveMd='extracted';


rv={'rawgaeCoh'};

behav={'exploring'};
% x axis limits
xl=[2 16];

% --- prepare graphics
labelscale('fontSz',6,'scaleFac',.25,'lineW',1.5,'markSz',8);
ornt='portrait';
figdir='d:\projects\hippocampus\rmouse\paper_atropine\rawFig\';
printas='-dpsc2';
% printas=[];

figName=mfilename;

figure(1), clf, orient(ornt);
subplot(1,3,1)

fi=1;
% --- paths, channels
rmouse_ini;
cd([WP.rootPath '\beta3_wtko\wt0001_04708']);
a001_r1;

AP_beta3_wtko;
rmouse_APcheck;
rawCh=rmouse_chan;
% 7 channels up to princ chan
chixi=AP.pcIdx-6:AP.pcIdx;

if isempty(strfind(DS.dpath,':')), DS.dpath=[WP.rootPath DS.dpath]; end

% --- load
switch saveMd
  case 'original'
    load([AP.resPath '\' AP.resFn],'r');
    F=r(1).rawgaeF;
    % --- all about time: axes
    xax=F;
    cix=find(xax>=xl(1) & xax<=xl(2));
    xax=xax(cix);

    % --- extract
    % determine index to results obtained with specified behavior in r
    % note that intersect will sort alphabetically, which we dont want)
    [nada,rix]=intersect({r(:).segmentType},behav);
    for rvi=1:length(rv)
      c=[];
      for bi=sort(rix)
        eval(['tmpc=r(bi).' rv{rvi} ';']);
        % structure of variables: freq|channel|channel
        tmpc=tmpc(:,chixi,chixi(1));
        % cut down
        tmpc=tmpc(cix,:,:);
      end
      eval([rv{rvi} '=tmpc;']);
      ylab='Coherence';
      save([figdir DS.abfFn '_tmpDat'],'xax','cix','F',rv{rvi},'ylab');
    end

  case 'extracted'
    load([figdir DS.abfFn '_tmpDat']);
end

% --- plot
for bi=1:length(behav)
  hold on
  
  rawgaeCoh=rawgaeCoh+repmat((1:size(rawgaeCoh,2))*-.9,size(rawgaeCoh,1),1);
  ph=plot(xax,rawgaeCoh,'k');
  axis on

  niceyax 
  set(gca,'xtick',[0:4:20]);

end

if ~isempty(printas),
  print(printas,[figdir figName]);
end




