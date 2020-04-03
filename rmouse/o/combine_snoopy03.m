function d=combine_snoopy03
% ** function d=combine_snoopy03
% does a PC analysis of key parameters, comparing wt with ko
% - assumes presence of auto.mat and cross.mat
% - assumes that ANAPR and DSET in these files are identical

global ANPAR DSET

printas=[];'-djpeg90';

qq={,'auto','cross'};
% choose behaviors to be dealt with (legal value of AP.segmentType)
behav={'immobile','exploring'};
nBehav=length(behav);

% obtain results variable(s) to collect and average/plot
auto_rv=set_rv('auto');
cross_rv=set_rv('cross');

% the parameters to make it into PC: 
sel_auto_rv={'rawPMnPeak','rawPMnPeakT','rawGaPEMn',...
             'thgaeCCPeakMn','thgaeCCPeakTMn'...
             'thCCPosPeakDecayMn'};
sel_cross_rv={'thCCPeakMn','thCCPeakTMn',...
               'gaCCPeakMn','gaCCPeakTMn',...
               'rawCohMnTh','rawCohMnGa',...
               'gaeCCPeakMn','gaeCCPeakTMn'};

% sel_auto_rv={'rawPMnPeak','rawGaPEMn','thgaeCCPeakMn'};
% sel_cross_rv=[];
 
sel_rv=cat(2,sel_auto_rv,sel_cross_rv);
nAutoSelPar=length(sel_auto_rv);
nCrossSelPar=length(sel_cross_rv);
nSelPar=nAutoSelPar+nCrossSelPar;


% electrode positions (first column in R.d) for which to grab results: 
% - 0 for the auto results
% - .6 for the cross results EXCEPT gaeCCPeakMn & T
dDistArr=[repmat(0,1,nAutoSelPar) repmat(.6,1,nCrossSelPar)];
% watch out, this gets both 'gaeCCPeakMn' and 'gaeCCPeakTMn'
tmpix=strmatch('gaeCCPeak',sel_rv);
if ~isempty(tmpix)
  dDistArr(tmpix)=.3;
end

[foo,auto_ix]=intersect(auto_rv,sel_auto_rv);
[foo,cross_ix]=intersect(cross_rv,sel_cross_rv);

% generate index to be used for concatenated R.d (and other fields)
combin_ix=[auto_ix cross_ix+length(auto_rv)];


rmouse_ini;

% -------- PART I: collection of data
% it is important to load the ANPAR and DSET that generated the data in the matfile
% - assume that ANPAR and DSET in cross.mat and auto.mat are identical
load([WP.rootPath '\beta3_wtko\' qq{1} '.mat'],'ANPAR','DSET');

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
    ndset(ri,1,i3)=length([DSET(ri,:,i3).nsRng]);
  end
end

% the variable containing selected parameters for just the principal channel:
% - observations (=animals) down the columns
% - one column per parameter
% - behaviors in slices
pcd=repmat(nan,[sum(ndset(1,1,:)), nSelPar, nBehav]);
% - tags for genotype
gttag=[];
for i3=1:n3
  gttag=[gttag; repmat(i3,ndset(1,1,i3),1)];
end

for ii=1:length(qq)
  load([WP.rootPath '\beta3_wtko\' qq{ii} '.mat'],'R','bix');
  eval([qq{ii} '_R=R; clear R;']);
end
% now we should have auto_R and cross_R in the workspace
% concatenate their .d fields and .indv.
dd=cat(2,auto_R.d,cross_R.d);
indv=cat(2,auto_R.indv,cross_R.indv);

nAnimal=ndset(1,1,i3);
cumsanimix=squeeze(cumsum(ndset(1,1,:)));

% loop over genotypes=slices of R.d
for i3=1:n3
  % index into pcd for animals (=observations)
  if i3==1
    anix=1:cumsanimix(1);
  else
    anix=cumsanimix(i3-1)+1:cumsanimix(i3);
  end
  % loop over behaviors=rows of dd
  for bi=1:length(bix)
    % loop over parameters=columns of dd
    for sel_rvi=1:length(combin_ix)
      tmpr=dd{bi,combin_ix(sel_rvi),i3};
      % identify electrodes of interest - unfortunately, because there may be missing
      % electrodes, this has to be done animal by animal
      individuals=indv{bi,combin_ix(sel_rvi),i3};
      uindv=unique(individuals);
      % on the occasion, a check
      if length(uindv)~=ndset(1,1,i3), error('ksdfk'); end
      for uix=1:length(uindv)
        locTmpr=tmpr(individuals==uindv(uix),:);
        [foo,finalix]=min(abs(locTmpr(:,1)+dDistArr(sel_rvi)));
        pcd(anix(uix),sel_rvi,bi)=locTmpr(finalix,2);
      end
    end
  end
end

global D ND
D=pcd(:,:,2);

[pcs,ND]=PCexplore('plotType',{'scatter'},'tag',gttag,'nPC',5,'normalize',1,'omitMd','cha');
return
plot(D(1:7,1),D(1:7,2),'bo')
hold on
plot(D(8:14,1),D(8:14,2),'ro')