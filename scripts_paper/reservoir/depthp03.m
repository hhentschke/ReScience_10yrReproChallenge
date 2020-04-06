% Generates depth profile plots of single results parameter produced by
% rmouse->combine4b. Also performs statistical analysis. There is a choice 
% between different types of 2- or 3-way ANOVA. Recording site is always one 
% (the first) factor; the other factors may be drug, behavior, genotype.
% 
% Current logic: one data file contains two factors, recording site and a
% factor the level of which is in the file name. Data from two or more
% files will be combined, thus permitting two sets of 2-way ANOVA, e.g.
% first rec depth x drug and then rec depth x behavior.

% TO DO:
% - check whether fitting works
% - logics of plot generation is designed for comparison of drug conditions only
%   -> adapt to new scheme of multiple comparison


% -------------------------------------
% ------- I. Prelims ------------------
% -------------------------------------
% datafile name
DFN={'combine4b_drug_exploring_gaeCCPeakMn_cross',...
    'combine4b_drug_immobile_gaeCCPeakMn_cross'};

DFN={'combine4b_drug_exploring_gaCCPeakTMn_cross',...
    'combine4b_drug_immobile_gaCCPeakTMn_cross'};

DFN={'combine4b_drug_exploring_gaCCPeakMn_cross',...
    'combine4b_drug_immobile_gaCCPeakMn_cross'};

DFN={'combine4b_drug_exploring_thCCPeakMn_cross',...
    'combine4b_drug_immobile_thCCPeakMn_cross'};

DFN={'combine4b_drug_exploring_thCCPeakTMn_cross',...
    'combine4b_drug_immobile_thCCPeakTMn_cross'};

DFN={'combine4b_drug_exploring_gaePMnPeak_auto',...
    'combine4b_drug_immobile_gaePMnPeak_auto'};

DFN={'combine4b_drug_exploring_gaeThNarrowPEMn_auto',...
    'combine4b_drug_immobile_gaeThNarrowPEMn_auto'};

DFN={'combine4b_drug_exploring_gaeCCPeakTMn_cross',...
    'combine4b_drug_immobile_gaeCCPeakTMn_cross'};

DFN={'combine4b_drug_exploring_gaeThPEMn_auto',...
    'combine4b_drug_immobile_gaeThPEMn_auto'};
  
DFN={'combine4b_drug_exploring_thCCPosPeakDecayMn_auto',...
    'combine4b_drug_immobile_thCCPosPeakDecayMn_auto'};

% --------- theta gammaEnv CC ------------  
DFN={'combine4b_drug_exploring_thgaeCCPeakMn_auto',...
    'combine4b_drug_immobile_thgaeCCPeakMn_auto'};

DFN={'combine4b_drug_exploring_thgaeCCPeakTMn_auto',...
    'combine4b_drug_immobile_thgaeCCPeakTMn_auto'};

DFN={'combine4b_drug_exploring_thgaeCCPeakPhaseMn_auto',...
    'combine4b_drug_immobile_thgaeCCPeakPhaseMn_auto'};

DFN={'combine4b_drug_exploring_thgaeCCZScore_auto',...
    'combine4b_drug_immobile_thgaeCCZScore_auto'};

DFN={'combine4b_drug_exploring_thgaeCCZTestP_auto',...
    'combine4b_drug_immobile_thgaeCCZTestP_auto'};

% --------- power measurements ------------  
DFN={'combine4b_drug_exploring_rawPMnPeak_auto',...
    'combine4b_drug_immobile_rawPMnPeak_auto'};

DFN={'combine4b_drug_exploring_gaePMnPeak_auto',...
    'combine4b_drug_immobile_gaePMnPeak_auto'};

DFN={'combine4b_drug_exploring_rawGaPEMn_auto',...
    'combine4b_drug_immobile_rawGaPEMn_auto'};

DFN={'combine4b_drug_exploring_rawGaNarrowPEMn_auto',...
    'combine4b_drug_immobile_rawGaNarrowPEMn_auto'};

DFN={'combine4b_drug_exploring_rawBePEMn_auto',...
    'combine4b_drug_immobile_rawBePEMn_auto'};

DFN={'combine4b_drug_exploring_rawDePEMn_auto',...
    'combine4b_drug_immobile_rawDePEMn_auto'};

DFN={'combine4b_drug_exploring_rawThPEMn_auto',...
    'combine4b_drug_immobile_rawThPEMn_auto'};

DFN={'combine4b_drug_exploring_rawThNarrowPEMn_auto',...
    'combine4b_drug_immobile_rawThNarrowPEMn_auto'};

% --------- the one currently tested ------------  
DFN={'combine4b_drug_exploring_thCCPosPeakDecayMn_auto',...
    'combine4b_drug_immobile_thCCPosPeakDecayMn_auto'};

% DFN={'combine4b_drug_exploring_thNegPeakCvIPI_auto',...
%     'combine4b_drug_immobile_thNegPeakCvIPI_auto'};

DFN={'exploring_thGaComod_auto',...
  'immobile_thGaComod_auto'};


% which kinda plot - symbols + lines (of fitted funcs) or bars only?
ptype='symb'; 
ptype='bar';

% switch data sets? (may look better on plots)
swid=strcmpi(ptype,'bar');
% re-generate fits and plot them along with data?
genFit=strcmpi(ptype,'symb');
% ctrl vs drug: colors & symbols
pset={[.6 .6 .6],[0 0 0];'s','o'};
if swid, pset=fliplr(pset); end
% print? 
printas='-dpsc2';
printas=[];

% appearance of plots
% error bars of bar plots look best with lineW 1
% small insets
labelscale('fontSz',7,'scaleFac',.2,'lineW',.5,'markSz',10); 
% bar plots
labelscale('fontSz',16,'scaleFac',.6,'lineW',.5,'markSz',10); 
% this one for lag plots
labelscale('fontSz',16,'scaleFac',.6,'lineW',.7,'markSz',12); 
% this one for lag plots
labelscale('fontSz',16,'scaleFac',.6,'lineW',.6,'markSz',9); 

% nov 05: scheme above is inconsistent. bar plots look OK if scalefac=.6 
% but reduction in corel draw is set to 30%. Even thus, they are larger than 
% the original plots of fig 3. Also, the scheme relies on lines not being scaled down

rmouse_ini;
ornt='landscape';
figdir='d:\projects\hippocampus\rmouse\paper_atropine\rawFig\';
figName=mfilename;

% --- containers for concatenations of data from files and the various 
% grouping variables 
ds_c=[];
indv_c=[];
groupTag_c=[];

% -------------------------------------
% ------- II. Load data ------------------
% -------------------------------------

for k=1:length(DFN)
  dfn=DFN{k}
  load([WP.rootPath '\beta3_wtko\export\' dfn '.mat']);
  % check for missing electrodes or animals within file
  if ~isequal(ds1(:,1),ds2(:,1)) || ~isequal(indv1,indv2)
    error(['missing values for factor' compareFactor]);
  end
  % --- concatenations
  % - for data within present file
  ds12=[ds1; ds2];
  indv12=[indv1; indv2];
  % rec depth | data set (i.e. ds1 vs ds2)
  groupTag=[ds12(:,1) [zeros(size(indv1)); ones(size(indv2))]];
  % - for comparisons among files
  ds_c=cat(1,ds_c,ds12);
  indv_c=cat(1,indv_c,indv12);
  % rec depth | data set | k (=file)
  groupTag_c=cat(1,groupTag_c,[groupTag k*ones(size(indv12))]);
end
% check for missing electrodes or animals between files
ds1=ds_c(groupTag_c(:,3)==1,1);
ds2=ds_c(groupTag_c(:,3)==2,1);
indv1=indv_c(groupTag_c(:,3)==1);
indv2=indv_c(groupTag_c(:,3)==2);

if ~isequal(ds1(:,1),ds2(:,1)) || ~isequal(indv1,indv2)
  error(['missing values for factor' compareFactor]);
end

% the factors (=names of group tags) to be compared within the entire data
% set (=the data in all files listed in DFN above)
groupTagName_c={'rec depth',compareFactor,'other'};
% the parameter investigated is in the file name(s), so extract it
tmp=strfind(dfn,'_');
rv=dfn(tmp(end-1)+1:tmp(end)-1);
% the levels of the factors compared: within each data file, levels of the
% first factor, recording depth, are in the first column of mds1 and mds2.
% The second factor's levels are listed in variable compStr. In the
% comparison of results between files we have to extract the levels from
% the file names
% compareFactorLevel{1}=num2str(mds1(:,1));
% compareFactorLevel{2}=compStr;

% delete vars not needed anymore, including means (mds1, mds2) - they will
% be re-computed below, which is mandatory when individual rec sites or
% animals shall be purged
clear ds1 ds2 mds1 mds2 ds12 indv1 indv2 indv12 compStr

% -------------------------------------
% --- III. Statistics & plots -------------
% -------------------------------------
% *** for alternative statistics or rearrangement of data for
% export & analysis in SPSS see depthp02 ***

% set up data such that the following comparisons can be made
% i) control vs. drug for both behaviors
% ii) immobile vs. exploring for both drug conditions
for g=1:size(groupTag_c,2)
  gtUVal{g}=unique(groupTag_c(:,g));
end

gtnIx=[3 2]
% for same functionality as depthp02 set k=1:1
for k=1:1
  circK=mod(k,2)+1;
  groupTagName=groupTagName_c([1 gtnIx(circK)]);
  disp(['********* factors: ' groupTagName{1} ', ' groupTagName{2}]);
  for g=1:2
    % fixix points to the constant factor in current comparison
    fixix=groupTag_c(:,gtnIx(k))==gtUVal{gtnIx(k)}(g);
    ix1=groupTag_c(:,gtnIx(circK))==gtUVal{gtnIx(circK)}(1);
    ix2=groupTag_c(:,gtnIx(circK))==gtUVal{gtnIx(circK)}(2);    

    ds1=ds_c(fixix & ix1,:);
    ds2=ds_c(fixix & ix2,:);
    ds12=[ds1; ds2];
    indv1=indv_c(fixix & ix1,:);
    indv2=indv_c(fixix & ix2,:);
    indv12=[indv1; indv2];
    
    meth=3;
    switch meth
      case 3
        % II. 2-way ANOVA with both factors as repeated measures
        % & post-hoc ttests with simplest corrections for multiple comparisons
        uIndv=unique(indv1);
        nIndv=length(uIndv);
        uRecSite=unique(ds1(:,1));
        nRecSite=length(uRecSite);

        % kill rec depth 0.6 for ANOVA because we need a balanced design
        ix=ds12(:,1)<.6;

        % *** exclude certain individuals ?
        % ix=ix&indv12(:,1)~=4;

        % for cross measurements comparing the auto components does not make
        % sense because they will be 1
        if ~isempty(strfind(dfn,'cross'))
          ix=ix&ds12(:,1)>.0;
        end
        ds12_new=ds12(ix,:);
        indv12_new=indv12(ix,:);
        groupTag_new=groupTag(ix,:);

        stats = rm_anova2(ds12_new(:,2),indv12_new,groupTag_new(:,1),groupTag_new(:,2),groupTagName)

        % now, if there's a drug OR interaction effect, investigate all electrodes
        p=[];
        if stats{strmatch(groupTagName{2},stats(:,1),'exact'),end}<=.05 ||...
            stats{strmatch(['rec depth x ' groupTagName{2}],stats(:,1),'exact'),end}<=.05
          for gg=1:nRecSite
            ix=ds1(:,1)==uRecSite(gg);
            [h,p(gg)]=ttest(ds1(ix,2),ds2(ix,2));
          end
          p
        end
    end

    % ---- fit (again)?
    if genFit
      [ft_,fo_,st_,ds1ix,ds2ix,ds12ix]=curveFit2rmousePar(ds1,ds2,ds12,rv);

      % --- fit and determine quality of fit for wt and ko
      set(fo_,'Startpoint',st_);
      ds1f=fit(ds1(ds1ix,1),ds1(ds1ix,2),ft_ ,fo_);
      ds2f=fit(ds2(ds2ix,1),ds2(ds2ix,2),ft_ ,fo_);
      ds12f=fit(ds12(ds12ix,1),ds12(ds12ix,2),ft_ ,fo_);

      % --- part 4: create curves representing the fits & plot
      % 1. all data pts & fit
      fitx=ds12(1,1):(ds12(end,1)-ds12(1,1))/200:ds12(end,1);
      % cfit objects will be 'fevaluated' automatically
      ds1fit=ds1f(fitx);
      ds2fit=ds2f(fitx);
    end

    % finally, (re-) compute means and std
    for gg=1:nRecSite
      ix=ds1(:,1)==uRecSite(gg);
      mds1(gg,1:3)=[uRecSite(gg) mean(ds1(ix,2)) std(ds1(ix,2))];
      mds2(gg,1:3)=[uRecSite(gg) mean(ds2(ix,2)) std(ds2(ix,2))];
    end
    
    % --- plot -----
    if swid
      mds=mds1;
      dsfit=ds1fit;
      mds1=mds2;
      ds1fit=ds2fit;
      mds2=mds;
      ds2fit=dsfit;
    end

    figure(g), clf, orient(ornt), hold on

    switch ptype
      case 'symb';
        ph=plot(mds1(:,2),mds1(:,1),pset{2,1});
        set(ph,'color',pset{1,1},'markerfacecolor',pset{1,1});
        ph=plot(ds1fit,fitx,'-');
        set(ph,'color',pset{1,1});

        ph=plot(mds2(:,2),mds2(:,1),pset{2,2});
        set(ph,'color',pset{1,2},'markerfacecolor',pset{1,2});
        ph=plot(ds2fit,fitx,'-');
        set(ph,'color',pset{1,2});
        x1=mds1(:,1);
        x2=mds2(:,1);

      case 'bar'
        if strcmp(WP.mver(1:2),'7.')
          bh=barh('v6',mds1(:,1),[mds1(:,2) mds2(:,2)],1.0,'grouped');
        else
          bh=barh(mds1(:,1),[mds1(:,2) mds2(:,2)],1.0,'grouped');
        end
        set(bh(1),'facecolor',pset{1,1});
        set(bh(2),'facecolor',pset{1,2});
        % these are the midpoints of the bars, subsequently to be used as x values
        % for error bars and/or fits
%         x1=rot90(mean(get(bh(1),'YData')));
%         x2=rot90(mean(get(bh(2),'YData')));
        x1=(mean(get(bh(1),'YData')))';
        x2=(mean(get(bh(2),'YData')))';
    end

    errorcross([mds1(:,2) x1],mds1(:,3)*[1 0],'color',pset{1,1});
    errorcross([mds2(:,2) x2],mds2(:,3)*[1 0],'color',pset{1,2});

    if g==1
      axis tight
      xl=get(gca,'xlim');
    end
    if ~isempty(strfind(dfn,'rawPPeakTMn'));
      set(gca,'xlim',[6 10],'xtick',[0:2:12]);
    elseif ~isempty(strfind(dfn,'thCCPeakTMn'));
      set(gca,'xlim',[-10 82],'xtick',[0:20:90]);
    elseif ~isempty(strfind(dfn,'gaCCPeakTMn'));
      set(gca,'xlim',[-2 12],'xtick',[-2:2:12]);
      % this is used for overlay of gamma with gammaEnv
      % set(gca,'xlim',[-2 14],'xtick',[]);
    elseif ~isempty(strfind(dfn,'thgaeCCPeakTMn'));
      set(gca,'xlim',[-25 60],'xtick',[-20:20:60]);
    elseif ~isempty(strfind(dfn,'gaeCCPeakTMn'));
      set(gca,'xlim',[-2 12],'xtick',[-2:2:12]);
    elseif ~isempty(strfind(dfn,'thCCPosPeakDecayMn'));
      set(gca,'xlim',[0 .7],'xtick',[0:.3:.6]);
    elseif ~isempty(strfind(dfn,'thNegPeakCvIPI'));
      set(gca,'xlim',[0 .37],'xtick',[0:.1:.4]);
    elseif ~isempty(strfind(dfn,'gaeThPEMn'));
      set(gca,'xlim',[0 .012],'xtick',[0:.005:.012]);
    elseif ~isempty(strfind(dfn,'thGaComod'));
      set(gca,'xlim',[0 .205],'xtick',[0:.1:.2]);
    else
      set(gca,'xlim',[max(xl(1),0) xl(2)+diff(xl)*.05]);
    end

    % ** override settings above
    %   xl=[-25 63];
    %   set(gca,'xlim',xl,'xtick',[-20:20:60]);

    xl=get(gca,'xlim');
    set(gca,'xaxisloc','top');
    % set(gca,'ylim',[-.05 .65],'ytick',[0:.1:.6]);
    set(gca,'ylim',[-.05 .65],'ytick',[]);
    if ~isempty(strfind(dfn,'PeakTMn')) && strcmpi(ptype,'symb') ;
      lh=line([0 0],get(gca,'ylim'),'linestyle','--','linewidth',get(ph,'linewidth')*.5,'color','k');
    end
    % xlabel(dfn);
    % ylabel('Recording depth (mm)');
    % rexy('ax',gca,'xfac',.4*.4,'yfac',.8*.4);
    rexy('ax',gca,'xfac',.4,'yfac',.8);

    % place markers for p-values in graph
    if ~isempty(p)
%       % p has to be flipped because it represents the electrodes in reverse
%       % order (compared to the plots)
%       p=fliplr(p);
      %     xl=get(gca,'xlim');
      %     xti=get(gca,'xtick');
      xoffs=xl(1)+diff(xl)*.92;
      % start with uncorrected values
      tmpp=p;
      tmpp(p>.05)=nan;
      tmpp(p<=.05)=1;
      ph=plot(xoffs,mds1(:,1).*tmpp','ko');
      set(ph,'markersize',13,'linewidth',.4);
      % now the Bonferronis
      if ~isempty(strfind(dfn,'cross')),
        divFac=nRecSite-1;
      else
        divFac=nRecSite;
      end
      tmpp=p;
      tmpp(p>.05/divFac)=nan;
      tmpp(p<=.05/divFac)=1;
      ph=plot(xoffs,mds1(:,1).*tmpp','k*');
      set(ph,'markersize',11,'linewidth',.4);
      % set(gca,'xlim',xl,'xtick',xti)
      set(gca,'xlim',xl);
    end

    groupTagName{2}
    
    if ~isempty(printas),
      print(printas,[figdir dfn]);
    end
  end
end
