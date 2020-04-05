function R=combine4
% combines individual rmouse data sets and computes grand averages
% needs ANPAR and DSET (=concatenation of individual and matching(!) AP and DS)
% difference to other combine funcs: 
% - *** compares groups of animals with unmatched data (e.g. different genotypes) 
%       => expects 3D DSET and ANPAR ***
% - optionally saves extracted variables in files cross.mat or auto.mat (so that this function, 
%   combine4.m, or combine4b.m, can bypass the lengthy process of combining data from all
%   experiments)
% Note also: if DSET and ANPAR are not truly 3D, everything up to saving data will work, 
% but not the plotting part. 

global ANPAR DSET

% choose electrode depth profile - physical or functional
dp='phys';'func';
% choose auto- or cross-channel results (cross not computed for theta, gamma env corr)

q='auto';
q='cross';

% interpolate missing channels?
ipol=0;

% loading/saving compiled data: if loadFlag is nonzero and cross.mat or auto.mat
% exist, these will be loaded up, otherwise the data will be compiled from scratch.
% If writeflag is nonzero, cross.mat or auto.mat will be generated/overwritten
loadflag=1;
writeflag=1;

% do the manova statistics?
doStats=0;

printas=[];'-djpeg90';

% choose behaviors to be dealt with (legal value of AP.segmentType)
behav={'immobile','exploring'};
% ..and corresponding colors for plots (currently not used)
bpstr={'m';'g'};
% genotypes
gt={'WT';'KO'};
% ..and corresponding symbols and colors for plots
gtpstr={'bo';'rs'};

% obtain results variable(s) to collect and average/plot
[rv,rvix]=set_rv(q);

close all;
% good if single figures are to be printed
labelscale('fontSz',12,'scaleFac',1.0,'lineW',1.5,'markSz',8); 
% multiple small figs
labelscale('fontSz',8,'scaleFac',.4,'lineW',.75,'markSz',5); 

rmouse_ini;

% -------- PART I: collection of data
% it is important to load the ANPAR and DSET that generated the data in the matfile
if loadflag & exist([WP.rootPath '\beta3_wtko\' q '.mat'],'file')
  load([WP.rootPath '\beta3_wtko\' q '.mat'],'ANPAR','DSET');
end
% check contents of DSET and ANPAR
[n1 n2 n3]=size(ANPAR);
[dn1 dn2 dn3]=size(DSET);
if ~isequal([n1 n2 n3],[dn1 dn2 dn3]), error('DSET and ANPAR must have equal dimension'); end
if n1>1
  warndlg('function has not been tested for n(rows)>1 - breaking');
  return
end

% the number of columns with nonempty fields of DSET tells us how many data sets there are in 
% reality for each condition (that is, genotype/concentration). There may be empty elements 
% because DSET and ANPAR are struct arrays that may have been automatically expanded
% during their generation
for i3=1:n3
  for ri=1:n1
    % almost any field of DSET will do the trick since they're all indispensable (and
    % mostly nonempty)
    ndsets(ri,1,i3)=length([DSET(ri,:,i3).nsRng]);
  end
end

% struct holding collected results: R
tmplt=cell(length(behav),length(rv),n3);
tmplt(:)={[]};
% all of the following fields are 3d cell arrays:
% - row=behavior (in order listed above)
% - col=parameter (in order listed above)
% - slice=genotype (in order listed above)
% - each element of the cell array contains this type of data:
R.d=tmplt;          % collected data: 2d arr, 1st col electrode pos, 2nd col value
R.indv=tmplt;       % code for individual animal/session (needed for ANOVA with repeated measures and may be handy otherwise)
R.ueix=tmplt;       % 1d cell array; for each of the ue, these are the indices into the corresponding R.d
R.ga=tmplt;         % 'grand average': 2d arr, holding  mean|std|N
% this one's a 2D cell array, containing statistics for comparison among genotypes
R.bstat=tmplt(:,:,1); 

loadSuccess=0;
if loadflag & exist([WP.rootPath '\beta3_wtko\' q '.mat'],'file')
  load([WP.rootPath '\beta3_wtko\' q '.mat'],'R','bix');
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
          % (control + drug, if any) the number of all LFP channes is invariant
          if ri==1
            % template (single data set)
            tempo=repmat(nan,nAllLFPCh,1);
            % template (whole experiment)
            tempo2=repmat(tempo,1,n1);
          end
          % load results var..
          if isempty(strfind(AP.resPath,':')), AP.resPath=[WP.rootPath AP.resPath]; end    
          load([AP.resPath '\' AP.resFn],'r');
          % ..find behaviors..
          for bi=1:length(behav)
            bix(bi)=strmatch(behav{bi},{r(:).segmentType});
          end
          % ..and extract
          for bi=1:length(bix)
            for rvi=1:length(rv)
              % tell what we're dealing with 
              disp(['data: ' AP.resFn ', behavior: ' behav{bi} ', data:' rv{rvi}]);      
              % tmpr is the original data 
              eval(['tmpr=r(bix(bi)).' rv{rvi} ';']);
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
                case {'thgaeCCPeakMn','thgaeCCPeakTMn','detheCCPeakMn','detheCCPeakTMn','thgaeCCPosPeakDecayMn','thgaeCCZScore','thgaeCCZTestP'}
                  if ~isempty(tmpr)
                    y=tempo;          
                    if strcmpi(rv{rvi},'thgaeCCZScore')
                      % Z>2.5 as criterion (corresponds to ~p=0.013): relative number of segments
                      nk=sum(isfinite(tmpr),1);
                      tmpr(~isfinite(tmpr))=0;
                      tmpr=sum(tmpr>2.5,1)./nk;
                    elseif strcmpi(rv{rvi},'thgaeCCZTestP')
                      % p<.05 as criterion: relative number of segments
                      nk=sum(isfinite(tmpr),1);
                      tmpr(~isfinite(tmpr))=1;
                      tmpr=sum(tmpr<.05,1)./nk;
                    end
                    y(AP.LFPccInd)=tmpr;
                  end
                otherwise
                  if length(diag(tmpr))>1
                    if strcmpi(q,'cross')
                      % extract CC data for principal channel: non analyzed channels=nan
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
                R.indv{bi,rvi,i3}(elix,1)=(i3-1)*100+ci;
              else
                % interpolate single missing electrodes??
                if ipol
                  % IMPORTANT: interpolation rests on the assumptions that
                  % - electrode spacing is equal among animals
                  % - depths are taken from a set of given values
                  % ==> electrode depths as specified in DS.rawCh must not have a gap
                  if strcmpi(dp,'func')
                    warndlg('interpolation of missing channel makes sense only if electrode depth profile (variable dp) is ''phys''');
                  else
                    % determine missing electrodes
                    if ~isempty(AP.LFPccOmitInd)
                      if length(AP.LFPccOmitInd)==1
                        disp(['interpolating missing channel:' DS.rawCh{AP.LFPccOmitInd,1}]);
                        figure(1), clf;
                        plot(WP.elx,y,'bo-');
                        hold on
                        y=interp1(WP.elx(AP.LFPccInd),y(AP.LFPccInd),WP.elx,'linear');
                        plot(WP.elx,y,'r+-');
                        drawnow;
                      else
                        warning([DS.abfFn ': more than one missing channel; no interpolation']);
                      end
                    end
                  end
                end % if:ipol
                % now insert
                R.d{bi,rvi,i3}(elix,ri+1)=y;
                % individual animal 'code': all wt=column order; all ko=100+column order 
                R.indv{bi,rvi,i3}(elix,1)=(i3-1)*100+ci;
              end % if:isempty(y)
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
        % find electrodes with any nan
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
      save auto
    case 'cross'
      save cross
  end
end


% ------ PART IV: statistics - two genotypes only, control only
if doStats
  for bi=1:length(bix)
    for rvi=1:length(rv)
      % union of all electrodes..
      uue=union(R.ue{bi,rvi,1},R.ue{bi,rvi,2});
      switch rv{rvi}
        case {'thgaeCCPeakMn','thgaeCCPeakTMn',...
            'detheCCPeakMn','detheCCPeakTMn',...
            'rawPMnPeak','rawPMnPeakT',...
            'rawDePEMn','rawThPEMn','rawGaPEMn','rawRiPEMn',...
            'thCCPosPeakDecayMn','thgaeCCPosPeakDecayMn'}
          % ..common electrodes at recording sites dorsal of & including SLM
          [cue,cueix1,cueix2]=intersect(R.ue{bi,rvi,1}(R.ue{bi,rvi,1}<=0), R.ue{bi,rvi,2}(R.ue{bi,rvi,2}<=0));
        otherwise
          % ..common electrodes at recording sites strictly dorsal of SLM
          [cue,cueix1,cueix2]=intersect(R.ue{bi,rvi,1}(R.ue{bi,rvi,1}<0), R.ue{bi,rvi,2}(R.ue{bi,rvi,2}<0));
      end
      % .. and indices to them in union of all
      [nix,supercix]=intersect(uue,cue);
      % matrix required as input to ANOVA-computing function rmaov2
      ad=[];
      % preallocate stats & set to nan
      % .P   p-value of H0 (no effect of genotype) from ANOVA
      % .p   1d array, p-value for each electrode
      R.bstat{bi,rvi}.p=nan*uue;
      % determine ue with data from EACH of the animals
      % do this by dealing with channel after channel because R.ueix contains cell arrays anyways
      ueix1=R.ueix{bi,rvi,1};
      ueix2=R.ueix{bi,rvi,2};
      for ei=1:length(cue)
        ix1=ueix1{cueix1(ei)};
        ix2=ueix2{cueix2(ei)};
        % bypass statistics for nonsense cases (principal electrodes where parameters are
        % all 0 or 1 (or any other common value) by definition (=due to normalization)
        if length(ix1)==ndsets(1,1,1) & length(ix2)==ndsets(1,1,2) % & length(unique([R.d{bi,rvi,1}(ix1,2);R.d{bi,rvi,2}(ix2,2)]))>1
          % genotype is first indep var, electrode depth second
          gt=1;
          ad=cat(1,ad,[R.d{bi,rvi,gt}(ix1,2)  repmat(gt,length(ix1),1)  R.d{bi,rvi,gt}(ix1,1)  R.indv{bi,rvi,gt}(ix1)]);
          gt=2;
          ad=cat(1,ad,[R.d{bi,rvi,gt}(ix2,2)  repmat(gt,length(ix2),1)  R.d{bi,rvi,gt}(ix2,1)  R.indv{bi,rvi,gt}(ix2)]);
          % R.bstat{bi,rvi}(supercix(ei))=p;
        end
      end
      disp(rv{rvi});
      disp(['electrodes: ' num2str((unique(ad(:,3)))')]);
      % **** 2-way ANOVA with repeated measures
      [R.bstat{bi,rvi}.P]=bwaov2(ad);
    end
  end
end

% % ------ PART IV: statistics - two genotypes only, control only
% for bi=1:length(bix)
%   for rvi=1:length(rv)
%     % union of all electrodes..
%     uue=union(R.ue{bi,rvi,1},R.ue{bi,rvi,2});
%     % .. common electrodes..
%     [cue,cueix1,cueix2]=intersect(R.ue{bi,rvi,1},R.ue{bi,rvi,2});    
%     % .. and indices to common electrodes in union of all
%     [nix,supercix]=intersect(uue,cue);
%     % preallocate stats & set to nan 
%     R.bstat{bi,rvi}=nan*uue;
%     % determine ue which have more than two values each
%     % do this by dealing with channel after channel because R.ueix contains cell arrays anyways
%     ueix1=R.ueix{bi,rvi,1};
%     ueix2=R.ueix{bi,rvi,2};    
%     for ei=1:length(cue)
%       ix1=ueix1{cueix1(ei)};
%       ix2=ueix2{cueix2(ei)};      
%       % bypass statistics for nonsense cases (principal electrodes where parameters are
%       % all 0 or 1 (or any other common value) by definition (=due to normalization) 
%       if length(ix1)>2 & length(ix2)>2 & length(unique([R.d{bi,rvi,1}(ix1,2);R.d{bi,rvi,2}(ix2,2)]))>2
%         [H,p]=ttest2(R.d{bi,rvi,1}(ix1,2),R.d{bi,rvi,2}(ix2,2),.05,'both');
%         % ranksum = Mann-Whitney U test => does not make sense with only a few data points
%         % [p,H]=ranksum(R.d{bi,rvi,1}(ix1,2),R.d{bi,rvi,2}(ix2,2),.05);        
%         R.bstat{bi,rvi}(supercix(ei))=p;             
%       end
%     end
%   end
% end

% ------ PART IV: plot (thus far, control condition only)
for rvi=1:length(rv)
  figure(rvi), orient tall
  for bi=1:length(bix)
    subplot(length(bix),1,bi), hold on
    title([behav{bi} ', ' rv{rvi} ', ' q]);
    uue=[];
    for i3=1:n3
      uue=union(uue,R.ue{bi,rvi,i3});      
      ph=errorbar(R.ue{bi,rvi,i3},R.ga{bi,rvi,i3}(:,1),zeros(size(R.ga{bi,rvi,i3}(:,1))),...
        R.ga{bi,rvi,i3}(:,2),gtpstr{i3});
    end
    % results from ANOVA:
    disp(rv{rvi});
    if R.bstat{bi,rvi}.P<.05
      disp(['***** H0 (identity of genotypes) rejected, p= ' num2str(R.bstat{bi,rvi}.P)]); 
      if R.bstat{bi,rvi}.P<.01, urtext('**',.9,'fontsize',30);
      else urtext('*',.9,'fontsize',30);
      end
    else
      disp(['H0 (identity of genotypes) not rejected, p= ' num2str(R.bstat{bi,rvi}.P)]);               
      urtext('n.s.');
    end
    
    % plot on top: crosses for nans, circles for p>=.05, single star for p<.05, double star for p<.01
    yl=get(gca,'ylim');
    yp=yl(2)-diff(yl)*[.05 .1];
    % nanix, refering to R.bstat, refers to uue also
    nanix=find(~isfinite(R.bstat{bi,rvi}.p));
    nsix=find(R.bstat{bi,rvi}.p>=.05);
    ix01=find(R.bstat{bi,rvi}.p<.01);  
    ix05=setdiff(find(R.bstat{bi,rvi}.p<.05),ix01);
    plot(uue(nanix), repmat(yp(1),size(nanix)),'kx');
    plot(uue(nsix), repmat(yp(1),size(nsix)),'ko'); 
    plot(uue(ix05), repmat(yp(1),size(ix05)),'k*');   
    plot(uue(ix01), repmat(yp(1),size(ix01)),'k*');     
    plot(uue(ix01), repmat(yp(2),size(ix01)),'k*');     
    grid on;
    niceyax;    
  end
  if ~isempty(printas), 
    print(printas,[WP.rootPath '\beta3_wtko\figures\' rv{rvi} '_' q '.jpg']); 
    % saveas(gcf,[WP.rootPath '\beta3_wtko\figures\' rv{rvi} '_' q ],'fig');
  end
end
    



