% generates depth profiles of all sorts of results vars
% warndlg('use depthp03??')
% datafile name

DFN={'combine4b_drug_exploring_gaeThPEMn_auto',...
    'combine4b_drug_immobile_gaeThPEMn_auto'};

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

DFN={'combine4b_drug_exploring_thCCPeakPhaseMn_cross_n4',...
    'combine4b_drug_immobile_thCCPeakPhaseMn_cross_n4'};

DFN={'combine4b_drug_exploring_thgaeCCPeakTMn_cross',...
    'combine4b_drug_immobile_thgaeCCPeakTMn_cross'};

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

DFN={'combine4b_drug_exploring_thNegPeakCvIPI_auto',...
    'combine4b_drug_immobile_thNegPeakCvIPI_auto'};

DFN={'combine4b_drug_immobile_thNegPeakCvA_auto',...
    'combine4b_drug_exploring_thNegPeakCvA_auto'};

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
DFN={'combine4b_drug_exploring_thCCPeakTMn_cross_n4',...
    'combine4b_drug_immobile_thCCPeakTMn_cross_n4'};

DFN={'combine4b_drug_exploring_thCCPeakPhaseMn_cross_n4',...
    'combine4b_drug_immobile_thCCPeakPhaseMn_cross_n4'};

DFN={'combine4b_drug_exploring_gaeCCPeakTMn_cross',...
    'combine4b_drug_immobile_gaeCCPeakTMn_cross'};

DFN={'combine4b_drug_exploring_gaeCCPeakPhaseMn_cross',...
    'combine4b_drug_immobile_gaeCCPeakPhaseMn_cross'};

DFN={'combine4b_drug_exploring_thgaeCCZScore_proportSignif_auto',...
    'combine4b_drug_immobile_thgaeCCZScore_proportSignif_auto'};

DFN={'combine4b_drug_exploring_thgaeCCZScore_auto',...
    'combine4b_drug_immobile_thgaeCCZScore_auto'};

DFN={'combine4b_drug_exploring_thgaeCCPeakMn_auto',...
    'combine4b_drug_immobile_thgaeCCPeakMn_auto'};

DFN={'combine4b_drug_exploring_rawThNarrowPEMn_auto',...
    'combine4b_drug_immobile_rawThNarrowPEMn_auto'};
% what kind of comparison - behavs or genotypes??
% *** the newer versions of data files exported by combine4b and loaded up
% here also contains variable compareFactor, overwriting its value set here
% ***
compareFactor='genotype';
compareFactor='behavior';
compareFactor='drug';

% which kinda plot - symbols + lines (of fitted funcs) or bars only?
ptype='symb'; 
ptype='bar';

% switch data sets? (may look better on plots)
swid=strcmpi(ptype,'bar');
% print? 
printas='-dpsc2';
printas=[];

% appearance of plots
% error bars of bar plots look best with lineW 1
% small insets
labelscale('fontSz',7,'scaleFac',.2,'lineW',.5,'markSz',10); 
% this one for lag plots
labelscale('fontSz',16,'scaleFac',.6,'lineW',.7,'markSz',12); 
% this one for lag plots
labelscale('fontSz',16,'scaleFac',.6,'lineW',.6,'markSz',9); 
% bar plots
labelscale('fontSz',16,'scaleFac',.6,'lineW',.5,'markSz',10); 

% nov 05: scheme above is inconsistent. bar plots look OK if scalefac=.6 
% but reduction in corel draw is set to 30%. Even thus, they are larger than 
% the original plots of fig 3. Also, the scheme relies on lines not being scaled down


rmouse_ini;
ornt='landscape';
figdir='d:\projects\rmouse\paper_atropine\rawFig\';
figName=mfilename;

switch compareFactor
  case 'genotype'
    % wt vs ko: colors & symbols
    pset={'b','r';'o','s'};
    % the pale versions
    % pset={[.85 .85 1],[1 .85 .85];'s','o'};
  case 'behavior'
    % exploring vs immobile: colors & symbols
    pset={[.6 .4 .1],[.5 1 .5];'s','o'};
  case 'drug'
    % ctrl vs drug: colors & symbols
    pset={[.6 .6 .6],[0 0 0];'s','o'};
end


% --- containers for concatenations of data from files and the various 
% grouping variables (in case we want to compare across files at the end
% of the loop)

ds12_c=[];
indv12_c=[];
groupTag_c=[];

groupTagName={'rec depth',compareFactor};
groupTagName_c={'rec depth',compareFactor,'other'};

% --- ** loop over files
for k=1:length(DFN)
  dfn=DFN{k}
  load([WP.rootPath '\beta3_wtko\export\' dfn '.mat']);

  
  % ***************************************************
  % ***************************************************
%   ds1(:,2)=abs(ds1(:,2));
%   ds2(:,2)=abs(ds2(:,2));  
% 
%   mds1(:,2)=abs(mds1(:,2));
%   mds2(:,2)=abs(mds2(:,2));  
  % ***************************************************
  % ***************************************************
  
  % --- concatenations
  % - for data within present file
  ds12=[ds1; ds2];
  indv12=[indv1; indv2];
  % rec depth | data set (i.e. ds1 vs ds2)
  groupTag=[ds12(:,1) [zeros(size(indv1)); ones(size(indv2))]];
  % - for comparisons among files
  ds12_c=cat(1,ds12_c,ds12);
  indv12_c=cat(1,indv12_c,indv12);
  % rec depth | data set | k
  groupTag_c=cat(1,groupTag_c,[groupTag k*ones(size(indv12))]);

  meth=0;
  switch meth
    case 0
      % forget all the restructuring above and arrange data in a way
      % needed for SPSS (2-way ANOVA with repeated measures).
      if ~isequal(ds1(:,1),ds2(:,1)) || ~isequal(indv1,indv2)
        error('we have a problem')
      end
      uIndv=unique(indv1);
      nIndv=length(uIndv);
      uRecSite=unique(ds1(:,1));
      nRecSite=length(uRecSite);
      % preallocation: as many rows as animals, number of columns =
      % nRecSite*2
      dsNew=repmat(nan,nIndv,2*nRecSite);
      for gg=1:nIndv
        iix=find(indv1==uIndv(gg));
        for hh=1:nRecSite
          rsix=find(ds1(:,1)==uRecSite(hh));
          superix=intersect(iix,rsix);
          if isempty(superix)
            disp(['no data found for animal #' int2str(uIndv(gg)) ', rec depth ' num2str(uRecSite(hh))]);
          else
            dsNew(gg,[hh hh+nRecSite])=[ds1(superix,2) ds2(superix,2)];
          end
        end
      end
      save(['d:\projects\rmouse\paper_atropine\spss\' dfn '.txt'],'dsNew','-ascii');

    case 1
      % 2-way within and between subjects ANOVA - doesn't work
      [p1,p2,p3]=bwaov2([ds12(:,2) groupTag(:,[1 2]) indv12]);
    case 2
      % ----- 2-way anova (electrode depth plus other comparison factor), no
      % repeated measures
      % define groups
      gT=cell(size(ds12,1),1);
      gT(1:size(ds1,1))=compStr(1);
      gT(size(ds1,1)+1:size(ds12,1))=compStr(2);
      gT={[ds1(:,1); ds2(:,1)],gT};

      [p,t,stats]=anovan([ds1(:,2);ds2(:,2)],gT,'display','off');
      p
      % c=multcompare(stats,'display','on','dimension',[1 2])
      % return
      
    case 3
      % II. 2-way ANOVA with both factors as repeated measures 
      % & post-hoc ttests with simplest corrections for multiple comparisons
      if ~isequal(ds1(:,1),ds2(:,1)) || ~isequal(indv1,indv2)
        error('we have a problem')
      end
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
      if stats{strmatch('drug',stats(:,1),'exact'),end}<=.05 ||...
        stats{strmatch('rec depth x drug',stats(:,1),'exact'),end}<=.05
        for gg=1:nRecSite
          ix=ds1(:,1)==uRecSite(gg);
          [h,p(gg)]=ttest(ds1(ix,2),ds2(ix,2));
        end
        p
      end
      
    case 4
      % ttest & friedman of paired samples - one rec site only
      % kill all rec depths except  0
      rd=0.1;
      ix=ds1(:,1)==rd;
      ds1=ds1(ix,:);
      ix=ds2(:,1)==rd;
      ds2=ds2(ix,:);
      [h,p1] = ttest(ds1(:,2),ds2(:,2));
      p2=friedman([ds1(:,2),ds2(:,2)],1,'off');
      disp(['ttest: ' num2str(p1) '; friedman: ' num2str(p2)]);
      
    case 5
      % average results from groups of electrodes representing e.g. the
      % dipoles (i.e., recording sites 0 and 0.1, and 0.5 and 0.6)
      avgRecSiteGroup={[0.0 0.1],[0.5 0.6]};
      % avgRecSiteGroup={[0.0 ],[0.3], [0.5 ]};
      nAvgRecSiteGroup=length(avgRecSiteGroup);
      % I. reshape data
      if ~isequal(ds1(:,1),ds2(:,1)) || ~isequal(indv1,indv2)
        error('we have a problem')
      end
      uIndv=unique(indv1);
      nIndv=length(uIndv);
      uRecSite=unique(ds1(:,1));
      nRecSite=length(uRecSite);
      % preallocation: as many rows as animals, number of columns =
      % nRecSite
      dsNew1=repmat(nan,nIndv,nRecSite);
      dsNew2=repmat(nan,nIndv,nRecSite);
      % this one will hold the averages: first rows control, then drug
      dsNew=repmat(nan,nIndv*2,nAvgRecSiteGroup);
      for gg=1:nIndv
        iix=find(indv1==uIndv(gg));
        for hh=1:nRecSite
          rsix=find(ds1(:,1)==uRecSite(hh));
          superix=intersect(iix,rsix);
          if isempty(superix)
            disp(['no data found for animal #' int2str(uIndv(gg)) ', rec depth ' num2str(uRecSite(hh))]);
          else
            dsNew1(gg,hh)=ds1(superix,2);
            dsNew2(gg,hh)=ds2(superix,2);
          end
        end
      end
      % now compute averages
      for gg=1:nAvgRecSiteGroup
        [nada,rsix]=intersect(uRecSite,avgRecSiteGroup{gg});
        dsNew(:,gg)=[nanmean(dsNew1(:,rsix),2); nanmean(dsNew2(:,rsix),2)];
      end
      
      % reshape again, into 1-column array, for input into rm_anova2
      dsNew=reshape(dsNew,nIndv*2*nAvgRecSiteGroup,1);
      indvTag=repmat([1:nIndv]',2*nAvgRecSiteGroup,1);
      % tag for recSite | drug 
      % ** watch out: tags for rec sites are 1, 2, ... whereas tags for
      % drug are 0 (control) and 1 **
      groupTag=makecol(ones(nIndv*2,1)*(1:nAvgRecSiteGroup));
      % groupTag=[zeros(nIndv*2,1); ones(nIndv*2,1)];
      groupTag(:,2)=repmat([zeros(nIndv,1); ones(nIndv,1)],nAvgRecSiteGroup,1);

      stats = rm_anova2(dsNew,indvTag,groupTag(:,1),groupTag(:,2),groupTagName)
      
      % finally, t-tests comparing drug effects ONLY within sites
      p=[];
      for gg=1:nAvgRecSiteGroup
        [h,p(gg)]=ttest(dsNew(groupTag(:,1)==gg & groupTag(:,2)==0),...
          dsNew(groupTag(:,1)==gg & groupTag(:,2)==1));
      end
      p
      
  end

  
  if swid
    if k==1
      pset=fliplr(pset);
    end
    mds=mds1;
    dsfit=ds1fit;
    mds1=mds2;
    ds1fit=ds2fit;
    mds2=mds;
    ds2fit=dsfit;
  end
  
  figure(k), clf, orient(ornt), hold on
  
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
      x1=rot90(mean(get(bh(1),'YData')));
      x2=rot90(mean(get(bh(2),'YData')));
  end
  
  errorcross([mds1(:,2) x1],mds1(:,3)*[1 0],'color',pset{1,1});
  errorcross([mds2(:,2) x2],mds2(:,3)*[1 0],'color',pset{1,2});
  
  if k==1
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
    set(gca,'xlim',[-25 63],'xtick',[-20:20:60]);
  elseif ~isempty(strfind(dfn,'thCCPeakTMn'));
    set(gca,'xlim',[-25 63],'xtick',[-20:20:60]);
  elseif ~isempty(strfind(dfn,'gaeCCPeakTMn'));
    % set(gca,'xlim',[-25 63],'xtick',[-20:20:60]);
    % insets
    set(gca,'xlim',[-2 12],'xtick',[-2:2:12]);
  elseif ~isempty(strfind(dfn,'thCCPosPeakDecayMn'));
    set(gca,'xlim',[0 .7],'xtick',[0:.3:.6]);
  elseif ~isempty(strfind(dfn,'thNegPeakCvA'));
    set(gca,'xlim',[0 1.05],'xtick',[0:.5:1]);
  elseif ~isempty(strfind(dfn,'thNegPeakCvIPI'));
    set(gca,'xlim',[0 .37],'xtick',[0:.1:.4]);
  elseif ~isempty(strfind(dfn,'gaeThPEMn'));
    set(gca,'xlim',[0 .012],'xtick',[0:.005:.012]);
  elseif ~isempty(strfind(dfn,'thGaComod'));
    set(gca,'xlim',[0 .205],'xtick',[0:.1:.2]);
  elseif ~isempty(strfind(dfn,'thgaeCCPeakPhaseMn')) || ~isempty(strfind(dfn,'thCCPeakPhaseMn'));
    set(gca,'xlim',[-1.25  3.5 ],'xtick',[-1:1:3]);
  elseif ~isempty(strfind(dfn,'thgaeCCZScore_proportSignif'));
    set(gca,'xlim',[0  1.2 ],'xtick',[0:.5:1]);
    
  else
    set(gca,'xlim',[max(xl(1),0) xl(2)+diff(xl)*.05]);  
  end
  
  


  xl=get(gca,'xlim');
  set(gca,'xaxisloc','top');
  % set(gca,'ylim',[-.05 .65],'ytick',[0:.1:.6]);
  set(gca,'ylim',[-.05 .65],'ytick',[]);
  if (~isempty(strfind(dfn,'PeakTMn')) || ~isempty(strfind(dfn,'PeakPhaseMn')))...
      && strcmpi(ptype,'symb') ;
    lh=line([0 0],get(gca,'ylim'),'linestyle','--','linewidth',get(ph,'linewidth')*.5,'color','k');
  end
  % xlabel(dfn);
  % ylabel('Recording depth (mm)');
  % rexy('ax',gca,'xfac',.4*.4,'yfac',.8*.4);
  rexy('ax',gca,'xfac',.4,'yfac',.8);
  
  % place markers for p-values in graph
  if ~isempty(p) 
    % p has to be flipped because it represents the electrodes in reverse
    % order (compared to mds1 or mds2)
    p=fliplr(p);
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
  
  if ~isempty(printas), 
    print(printas,[figdir dfn]); 
  end
end

