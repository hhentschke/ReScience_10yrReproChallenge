function R=combine4b
% this is a q&d modification of combine4b
% the specific job here: re-compute thgaeCCPeakMn and thgaeCCPeakTMn from
% segments with significant peak CC

% combines individual rmouse data sets and computes grand averages
% needs ANPAR and DSET (=concatenation of individual and matching(!) AP and DS)
% difference to other combine funcs: 
% - *** compares groups of animals with unmatched data (e.g. different genotypes) 
%       => expects 3D DSET and ANPAR ***
% - a few internal variables restructured/renamed
% difference to combine4: 
% - different statistics (F-test, Motulsky, necessitating fitting functions to data)



% improvements:

global ANPAR DSET

% which factor to compare?
% - for each type of behavior, genotypes will be compared
compareFactor='genotype'; 
% - within each genotype, behavioral episodes will be compared with each other
compareFactor='behavior';
% - for each type of behavior, control will be compared with drug
compareFactor='drug'; 

% the project directory
projSubDir='\WTb3N265M';
projSubDir='\beta3_wtko';

% shall data & fits be exported? (one file per analysis parameter, for plot routines) 
export=1;
% print?
printas='-djpeg90';
printas=[];

curFigPath=[projSubDir '\figures'];

% choose electrode depth profile - physical or functional
dp='phys';'func';
% choose auto- or cross-channel results (cross not computed for theta, gamma env corr)
q='cross';
q='auto';

% loading/saving compiled data: if loadFlag is nonzero and cross.mat or auto.mat
% exist, these will be loaded up, otherwise the data will be compiled from scratch.
% If writeflag is nonzero, cross.mat or auto.mat will be generated/overwritten
loadflag=0;
writeflag=0;

% choose behaviors to be compared/plotted (legal value of AP.segmentType)
behav={'immobile','exploring'};
% ..and corresponding symbols and colors for plots (colors first line)
bpset={[.6 .4 .1],[.5 1 .5];'s','o'};
% genotypes
gt={'WT';'KO'};
% ..and corresponding symbols and colors for plots (colors first line)
gtpset={'b','r';'o','s'};
% drug conditions
dr={'control';'atropine'};
% ..and corresponding symbols and colors for plots (colors first line)
drpset={'b','m';'o','s'};
% mm, limits of electrode depth (inclusive; slm=0, dorsal ones negative, ventral ones positive)
% ** june 22: for some strange reason, electrode depths are not equally spaced
% (specifically, what should be exactly 0.7 is in fact a number slightly
% higher). Hence the funny depthLim.
depthLim=[-.70001 0];
depthLim=[-.6 0];

% obtain results variable(s) to collect and average/plot
[rv,rvix]=set_rv(q);

rv={'thgaeCCPeakTMn_significant'};
rvi=1;

close all;
% multiple small figs
labelscale('fontSz',8,'scaleFac',.4,'lineW',.75,'markSz',5); 
% good if single figures are to be printed
labelscale('fontSz',12,'scaleFac',1.0,'lineW',1.5,'markSz',8); 

rmouse_ini;

% -------- PART I: collection of data
% it is important to load the ANPAR and DSET that generated the data in the matfile
if loadflag & exist([WP.rootPath projSubDir '\' q '.mat'],'file')
  load([WP.rootPath projSubDir '\' q '.mat'],'ANPAR','DSET');
end
[n1 n2 n3]=size(ANPAR);
[dn1 dn2 dn3]=size(DSET);
if ~isequal([n1 n2 n3],[dn1 dn2 dn3]), error('DSET and ANPAR must have equal dimension'); end

% struct holding collected results: R
tmplt=cell(length(behav),length(rv),n3);
tmplt(:)={[]};
% all of the following fields are 3d cell arrays:
% - row=behavior (in order listed above)
% - col=parameter (in order listed above)
% - slice=genotype (in order listed above)
% - each element of the cell array contains this type of data:
R.d=tmplt;          % collected data: 2d arr, 1st col electrode pos, 2nd col value (drug exprmnts: 3rd+ cols = values ar var. concentrations)
R.indv=tmplt;       % code for individual animal/session (needed for ANOVA with repeated measures and may be handy otherwise)
R.ueix=tmplt;       % 1d cell array; for each of the ue, these are the indices into the corresponding R.d
R.ga=tmplt;         % 'grand average': 2d arr, holding  mean|std|N
% this one's a 2D cell array, containing statistics for comparison among genotypes
R.bstat=tmplt(:,:,1); % 1d array, for each electrode holding p-values etc 

loadSuccess=0;
if loadflag & exist([WP.rootPath projSubDir '\' q '.mat'],'file')
  load([WP.rootPath projSubDir '\' q '.mat'],'R','bix');
  if exist('R','var')
    loadSuccess=1;
  end
else
  % loop over data sets: 
  % one experiment per column, concentration down the columns, different types of animals
  % (e.g. genotypes) in different slices
  for i3=1:n3
    for ci=1:n2
      for ri=1:n1
        AP=ANPAR(ri,ci,i3);
        DS=DSET(ri,ci,i3);
        % if one of the (indispensable) fields of DS isempty, the corresponding struct
        % element does not correspond to a data set 
        if ~isempty(DS.rawCh)
          rawCh=rmouse_chan;
          % let's make the reasonable assumption that for all data won in one experiment
          % the number of all LFP channes is invariant
          if ri==1
            % template (single data set)
            tempo=repmat(nan,nAllLFPCh,1);
            % template (whole experiment)
            tempo2=repmat(tempo,1,n1);
          end
          
%           % if dpath does not contain a drive letter, pre-pend WP.rootPath
%           if isempty(strfind(DS.dpath,':')), DS.dpath=[WP.rootPath DS.dpath]; end
%           % extract si from abf file - if matfile exists, pick it instead of abf file
%           if exist([DS.dpath '\' DS.abfFn '.mat'],'file')
%             load([DS.dpath '\' DS.abfFn '.mat'],'abfi');
%           elseif exist([DS.dpath '\' DS.abfFn '.abf'],'file')
%             abfi=abfinfo([DS.dpath '\' DS.abfFn '.abf'],'verbose',0);
%           else
%             error([DS.dpath '\' DS.abfFn ' does not exist'])
%           end
%           % sampling interval
%           si=abfi.si;
          
          % load results var..
          if isempty(strfind(AP.resPath,':')), AP.resPath=[WP.rootPath AP.resPath]; end    
          load([AP.resPath '\' AP.resFn],'r');
          % ..find behaviors..
          tmpb=cell(1,length(r));
          [tmpb{:}]=deal(r.segmentType);
          for bi=1:length(behav)
            bix(bi)=strmatch(behav{bi},tmpb);
          end
          % ..and extract
          for bi=1:length(bix)
            for rvi=1:length(rv)
              % tell what we're dealing with 
              disp(['data: ' AP.resFn ', behavior: ' behav{bi} ', data:' rv{rvi}]);      
               
              
              
              if bi==1 & rvi==1 & ri==1
                % index into R.d rows
                elix=[1:nAllLFPCh]+size(R.d{bi,rvi,i3},1);
              end
              % first data set: set up electrode position in first column
              if ri==1
                if strcmpi(dp,'phys')          
                  R.d{bi,rvi,i3}(elix,1)=WP.elx;
                elseif strcmpi(dp,'func')
                  R.d{bi,rvi,i3}(elix,1)=WP.felx;            
                end
              end
              y=[];
              % select extraction method depending on results var chosen 
              switch rv{rvi}
                case {'thgaeCCPeakTMn_significant'}
                  % lag
                  tmprT=r(bix(bi)).thgaeCCPeakT;
                  % peak amplitude
                  % tmprA=r(bix(bi)).thgaeCCPeak;
                  if ~isempty(tmprT)
                    y=tempo;
                    % do it
                    
                    for g=1:size(tmprT,2)
                      % index to all segments with significant cc
                      ix=r(bix(bi)).thgaeCCZScore(:,g)>2.0;
                      if isempty(ix)
                        disp('z score limit lowered to 1.5');
                        ix=r(bix(bi)).thgaeCCZScore(:,g)>1.5;
                      end
                      tmpr(g)=mean(tmprT(ix,g));
                    end
                    y(AP.LFPccInd)=tmpr;
                  end

                
                case {'thgaeCCPeakMn','thgaeCCPeakTMn','detheCCPeakMn','detheCCPeakTMn','thgaeCCPosPeakDecayMn','thgaeCCZScore','thgaeCCZTestP'}
                  if ~isempty(tmpr)
                    y=tempo;          
                    if strcmpi(rv{rvi},'thgaeCCZScore')
                      % Z>2.5 as criterion (corresponds to ~p=0.013): relative number of segments
                      nk=size(tmpr,1);
                      tmpr=sum(tmpr>2.5,1)./nk;
                    elseif strcmpi(rv{rvi},'thgaeCCZTestP')
                      % p<.05 as criterion: relative number of segments
                      nk=size(tmpr,1);
                      tmpr=sum(tmpr<.05,1)./nk;
                    end
                    y(AP.LFPccInd)=tmpr;
                  end
                otherwise
                  if length(diag(tmpr))>1
                    if strcmpi(q,'cross')
                      % extract CC data for principal channel: non analyzed channels 
                      % (nans) will be ignored on the plot
                      y=cat(1, tmpr{1:AP.LFPpcInd2,AP.LFPpcInd2});
                      y=[y; cat(1, tmpr{AP.LFPpcInd2,AP.LFPpcInd2+1:end})];
                    elseif strcmpi(q,'auto')
                      y=cat(1, tmpr{AP.dixie});              
                    end
                  end
              end % switch
              % concatenate all data sets: 1st col electrode pos, 2nd+ col values for
              % concentration in usual order (currently: control, drug, recovery)
              if isempty(y)
                warning(['data to be extracted do not exist']);
                R.d{bi,rvi,i3}(elix,ri+1)=tempo;          
              else
                R.d{bi,rvi,i3}(elix,ri+1)=y;
              end % if:isempty(y)
              % individual animal 'code': all wt=column order; all ko=100+column order
              R.indv{bi,rvi,i3}(elix,1)=(i3-1)*100+ci;
            end % for:par
          end % for:behav
        end % if ~isempty(DS.rawCh)
      end % for:rows of ANPAR=concs
    end % for:cols of ANPAR=experiments
  end % for:slices of ANPAR=genotypes
end
% (r not needed anymore)
clear r

if ~loadSuccess
  % ------ PART II: find ue and get indices
  % (this is done separately in order to have grand averaging and statistics run independently of each other)
  for i3=1:n3
    for bi=1:length(bix)
      for rvi=1:length(rv)
        tmpr=R.d{bi,rvi,i3};
        % kick out all electrodes with any nan
        badix=find(any(~isfinite(tmpr),2));
        if ~isempty(badix),
          disp([int2str(length(badix)) ' entries are NaNs']);
          tmpr(badix,:)=[];
          R.d{bi,rvi,i3}(badix,:)=[];
          R.indv{bi,rvi,i3}(badix,:)=[];
        end
        % electrode positions
        ue=unique(tmpr(:,1));
        % which we need to keep (although within behaviors and analysis parameters
        % of one group of animals no variance is expected)
        R.ue{bi,rvi,i3}=ue;
        % preallocate R.ga
        % dimensions: #electrodes | #parameters (mean,std,N)  | #conc (control,drug,recovery)
        R.ga{bi,rvi,i3}=repmat(nan,[length(ue),3,n1]);
        % preallocation of R.bstat makes no sense here
        ix={};
        for ei=1:length(ue)
          ix{ei}=find(tmpr(:,1)==ue(ei));
        end
        R.ueix{bi,rvi,i3}=ix;
      end
    end
  end

  % ------ PART III: grand averages
  for i3=1:n3
    for bi=1:length(bix)
      for rvi=1:length(rv)
        % compute averages and std (nans had been kicked out before)
        tmpr=R.d{bi,rvi,i3};
        % the order of elements in R.ueix{bi,rvi,i3} corresponds to that in R.ue{bi,rvi,i3}
        for ei=1:length(R.ueix{bi,rvi,i3})
          ix=R.ueix{bi,rvi,i3}{ei};
          % columns: mean|std|N; slices: control|drug|recovery
          R.ga{bi,rvi,i3}(ei,1,:)=reshape(mean(tmpr(ix,2:end),1),[1 1 n1]);
          R.ga{bi,rvi,i3}(ei,2,:)=reshape(std(tmpr(ix,2:end),0,1),[1 1 n1]);
          R.ga{bi,rvi,i3}(ei,3,:)=repmat(length(ix),[1 1 n1]);
        end
      end
    end
  end
end % if ~loadSuccess

% save results up to here?
if writeflag
  switch q
    case 'auto'
      save auto R ANPAR DSET behav bix rv q
    case 'cross'
      save cross R ANPAR DSET behav bix rv q
  end
end

switch compareFactor
  case 'genotype'
    loopP=bix;
    pset=gtpset;
  case 'behavior'
    % there are only two genotypes
    loopP=1:length(gt);
    pset=bpset;    
  case 'drug'
    % 'drug' is very much like 'genotype'
    loopP=bix;
    pset=drpset;    
  otherwise
    error('illegal choice of ''compareFactor''');
end

if ischar(rvix), rvix=length(rv):-1:1; end
% ------ PART IV: curve fitting for the sake of statistics 
for rvi=rvix
  figure(rvi), orient landscape
  % loop either over behaviors (comparison of genotypes) or genotypes (comparison of behaviors)
  for loopi=1:length(loopP)  
    % --- part 1: collect data, restrict to recording sites dorsal of & including SLM and
    % invert x axis so that independent var has values >=0
    % - data sets:
    switch compareFactor
      case 'genotype'
        ds1=R.d{loopi,rvi,1};
        ds2=R.d{loopi,rvi,2};
      case 'behavior'        
        ds1=R.d{1,rvi,loopi};
        ds2=R.d{2,rvi,loopi};
      case 'drug'        
        ds1=R.d{loopi,rvi,1}(:,[1 2]);
        ds2=R.d{loopi,rvi,1}(:,[1 3]);
    end
    % restrict depth range
    ds1(ds1(:,1)<depthLim(1) | ds1(:,1)>depthLim(2),:)=[];
    ds2(ds2(:,1)<depthLim(1) | ds2(:,1)>depthLim(2),:)=[];
    % invert sign
    ds1(:,1)=ds1(:,1)*-1;
    ds2(:,1)=ds2(:,1)*-1;
    % sort
    ds1=sortrows(ds1,1);
    ds2=sortrows(ds2,1);
    % - combination of both
    ds12=sortrows([ds1; ds2]);
    
    % for plots of average data points (that is, the original grand averages) do the same
    % thing (+shuffle & combine columns so that the resulting var has pos|mean|std cols):
    % data sets:
    switch compareFactor
      case 'genotype'
        mds1=[R.ue{loopi,rvi,1} R.ga{loopi,rvi,1}(:,1:2)];
        mds2=[R.ue{loopi,rvi,2} R.ga{loopi,rvi,2}(:,1:2)];
        titl=[behav{loopi} ', ' rv{rvi} ', ' q];
        fnp=[behav{loopi} '_' rv{rvi} '_' q];
        compStr=gt;
      case 'behavior'        
        mds1=[R.ue{1,rvi,loopi} R.ga{1,rvi,loopi}(:,1:2)];
        mds2=[R.ue{2,rvi,loopi} R.ga{2,rvi,loopi}(:,1:2)];
        titl=[gt{loopi} ', ' rv{rvi} ', ' q];
        fnp=[gt{loopi} '_' rv{rvi} '_' q];
        compStr=behav;
      case 'drug'
        mds1=[R.ue{loopi,rvi,1} R.ga{loopi,rvi,1}(:,1:2,1)];
        mds2=[R.ue{loopi,rvi,1} R.ga{loopi,rvi,1}(:,1:2,2)];
        titl=['wt, comprsn drug, ' behav{loopi} ', ' rv{rvi} ', ' q];
        fnp=[behav{loopi} '_' rv{rvi} '_' q]
        compStr=dr;
    end
    % restrict depth range
    mds1(mds1(:,1)<depthLim(1) | mds1(:,1)>depthLim(2),:)=[];
    mds2(mds2(:,1)<depthLim(1) | mds2(:,1)>depthLim(2),:)=[];
    % invert sign
    mds1(:,1)=mds1(:,1)*-1;
    mds2(:,1)=mds2(:,1)*-1;
    
    % --- part 2: depending on parameter transform data and set up model
    [ft_,fo_,st_,ds1ix,ds2ix,ds12ix]=curveFit2rmousePar(ds1,ds2,ds12,rv{rvi});
    
    % --- part 3: fit and determine quality of fit for wt and ko
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
    ds12fit=ds12f(fitx);
    % -----------------interlude: export data?
    if export
      fn=[mfilename '_' compareFactor(1:4) '_' fnp];
      save([WP.rootPath projSubDir '\export\' fn],'compareFactor','compStr','ds1','ds2','ds12','mds*','ds1fit','ds2fit','fitx');
    end
    % -------------------------------------
    subplot(length(loopP),2,(loopi-1)*2+1), hold on
    title(titl);
    % first
    ph=plot(ds1(:,1),ds1(:,2),pset{2,1});
    set(ph,'color',pset{1,1});
    ph=plot(fitx,ds1fit,'-');
    set(ph,'color',pset{1,1});    
    % second
    ph=plot(ds2(:,1),ds2(:,2),pset{2,2});
    set(ph,'color',pset{1,2});
    ph=plot(fitx,ds2fit,'-');
    set(ph,'color',pset{1,2});        
    % combo fit in black
    plot(fitx,ds12fit,'k-');
    nicexyax;
    % 2. averages & fit
    subplot(length(loopP),2,(loopi-1)*2+2), hold on
    ph=errorbar(mds1(:,1),mds1(:,2),mds1(:,3),pset{2,1});
    set(ph,'color',pset{1,1});
    ph=plot(fitx,ds1fit,'-');
    set(ph,'color',pset{1,1});    
    ph=errorbar(mds2(:,1),mds2(:,2),mds2(:,3),pset{2,2});
    set(ph,'color',pset{1,2});
    ph=plot(fitx,ds2fit,'-');
    set(ph,'color',pset{1,2});    
    niceyax;
    % --- part 5: F-test
    % statistical test for similarity (order: wt crf, ko crf, combo crf)
    [p,F,radj1,radj2]=curvecomp([ds1(ds1ix,:) ds1f(ds1(ds1ix,1))], ...
                                [ds2(ds2ix,:) ds2f(ds2(ds2ix,1))],...
                                [ds12(ds12ix,:) ds12f(ds12(ds12ix,1))],...
                                length(ds1ix)-length(st_),length(ds2ix)-length(st_));

    urtext(['p=' num2str(p,'%1.3f')],.85,'fontsize',12);
    if p<.05
      disp(['***** H0 (identity of depth-response profiles) rejected, p= ' num2str(p)]); 
%       if p<.01, urtext('**',.9,'fontsize',20);
%       else urtext('*',.9,'fontsize',20);
%       end
    else
      disp(['H0 (identity of depth-response profiles) not rejected, p= ' num2str(p)]);               
%       urtext('n.s.');
    end
    % information on goodness of fit
    ultext(['R_{adj}(1):' num2str(radj1,3) '; R_{adj}(2):' num2str(radj2,3)]);
    ds1f
    ds2f
    ds12f
  end % for: behaviors
  if ~isempty(printas)
    print(printas,[WP.rootPath projSubDir '\figures\' compareFactor(1:4) '_' rv{rvi} '_' q '_fit.jpg']); 
    % saveas(gcf,[WP.rootPath projSubDir '\figures\' rv{rvi} '_' q '_fit'],'fig');
  end
end % for: parameters





        
%       case {'thgaeCCPeakTMn_OK'}
%         % including data for x=0 does make sense for this particular CC!
%         ds1ix=find(ds1(:,1)>=0);
%         ds2ix=find(ds2(:,1)>=0);
%         ds12ix=find(ds12(:,1)>=0);
%         % Hill with open upper limit, offset and time shift; a=Hill coeff, b=EC50, m=max, o=offset,
%         % d=time offset
%         ft_ = fittype('m*(x+d)^a/(b^a + (x+d)^a)+o' ,...
%           'dependent',{'y'},'independent',{'x'},...
%           'coefficients',{'a','b','m','o','d'});
%         fo_ = fitoptions('method','NonlinearLeastSquares','Lower',[0 .1 20 -50 .001],'Upper',[10000 .8 100 30 1]);
%         % starting values for parameters - potentially sensitive!
%         st_ = [1 .4 40 -10 .01];
%         
%       case {'thgaeCCPeakTMn_OK2'}
%         % including data for x=0 does make sense for this particular CC!
%         ds1ix=find(ds1(:,1)>=0);
%         ds2ix=find(ds2(:,1)>=0);
%         ds12ix=find(ds12(:,1)>=0);
%         % Hill with open upper limit and offset; a=Hill coeff, b=EC50, m=max, o=offset
%         ft_ = fittype('m*x^a/(b^a + x^a)+o' ,...
%           'dependent',{'y'},'independent',{'x'},...
%           'coefficients',{'a','b','m','o'});
%         fo_ = fitoptions('method','NonlinearLeastSquares','Lower',[0 .1 20 -50],'Upper',[10000 .8 100 30]);
%         % starting values for parameters - potentially sensitive!
%         st_ = [1 .4 40 -10];
% 
%       case {'thCCPeakTMn_OK'}
%         % including data for x=0 does not make sense for CC parameters
%         ds1ix=find(ds1(:,1)>0);
%         ds2ix=find(ds2(:,1)>0);
%         ds12ix=find(ds12(:,1)>0);
%         % Hill with upper limit open; a=Hill coeff, b=EC50, m=max
%         ft_ = fittype('m*x^a/(b^a + x^a)' ,...
%           'dependent',{'y'},'independent',{'x'},...
%           'coefficients',{'a', 'b', 'm'});
%         fo_ = fitoptions('method','NonlinearLeastSquares','Lower',[0 .1 20],'Upper',[10000 .8 100]);
%         % starting values for parameters - potentially sensitive!
%         st_ = [1 .4 60];
%         

%       case {'thCCPeakMn'}
%         % including data for x=0 does not make sense for CC parameters
%         ds1ix=find(ds1(:,1)>0);
%         ds2ix=find(ds2(:,1)>0);
%         ds12ix=find(ds12(:,1)>0);
%         % sum of two gaussians + offset, force first hump to reside at x=0 with amplitude of 1
%         % and second hump + offset to cancel each other at x=0
%         % -> nice try, fit looks generally OK; serious theoretical drawback is that the
%         % larger hump of the composite is always at x>0 and has amplitude >1
%         ft_ = fittype('exp(-b1*(x)^2) + a2*exp(-b2*(x-c2)^2) - a2*exp(-b2*(c2)^2)' ,...
%           'dependent',{'y'},'independent',{'x'},...
%           'coefficients',{'b1','a2','b2','c2'});
%         fo_ = fitoptions('method','NonlinearLeastSquares');
%         % fo_ = fitoptions('method','NonlinearLeastSquares',...
%           % 'Lower',[.0001 .0001 -50 -50],'Upper',[1000 100 50 100]);
%         % starting values for parameters
%         st_ = [10 .6 5 .7];




%       case {'gaeCCPeakTMn'}
%         % including data for x=0 does not make sense for CC parameters
%         ds1ix=find(ds1(:,1)>0);
%         ds2ix=find(ds2(:,1)>0);
%         ds12ix=find(ds12(:,1)>0);
%         % gaussian * sin (crap)
%         ft_ = fittype('exp(-b1*(x-c1)^2) * a2*sin(b2*x)' ,...
%           'dependent',{'y'},'independent',{'x'},...
%           'coefficients',{'b1','c1','a2','b2'});
%         fo_ = fitoptions('method','NonlinearLeastSquares');
%         % fo_ = fitoptions('method','NonlinearLeastSquares',...
%           % 'Lower',[.0001 .0001 -50 -50],'Upper',[1000 100 50 100]);
%         % starting values for parameters - it is of utmost importance that b2, the cosine
%         % freq, has a good starting value
%         st_ = [20 .25 .25 6];

%       case {'gaeCCPeakTMn'}
%         % including data for x=0 does not make sense for CC parameters
%         ds1ix=find(ds1(:,1)>0);
%         ds2ix=find(ds2(:,1)>0);
%         ds12ix=find(ds12(:,1)>0);
%         % not too bad, but the periodicity of the sine looks a little funnny at the edges
%         ft_ = fittype('a1*exp(a2*(sin(c2*x))^4) + a3*x - a1' ,...
%           'dependent',{'y'},'independent',{'x'},...
%           'coefficients',{'a1','a2','c2','a3'});
%         fo_ = fitoptions('method','NonlinearLeastSquares');
%         % fo_ = fitoptions('method','NonlinearLeastSquares',...
%           % 'Lower',[.0001 .0001 -50 -50],'Upper',[1000 100 50 100]);
%         % starting values for parameters
%         st_ = [.1 3 5 -1];
%         
