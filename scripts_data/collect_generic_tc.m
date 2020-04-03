% collect DSETs and ANPARs that belong to a single atropine session 
% ** NOTE:
% - script collects AP with the 'standard' segment length of ~4 s
% - purpose of the exercise is to generate overview plots with combine_tc
%   or combine_bv; since neither of these do any statistics we don't need the
%   'compareFactor' variable here

global ANPAR DSET WP

% ANPAR and DSET will be set up as 2D struct arrays. First row is control, second+
% rows are post-drug-administration files, one column per experiment
ANPAR=[];
DSET=[];

% learn data root path from rmouse ini routine
rmouse_ini;
ci=0;

% select experiments
% two options:
% (i)  select single experiment, then generate time course plot of selected parameters
%      by manually calling combine_tc (which accepts only single columns)
% (ii) select all experiments with behavioral scoring, time course of behaviors will be 
%      automatically generated via combine_bv (see below)

expr={'wt2642','wt2654','wt2008','wt0002','wt0001','wt0003'};
expr={'wt2642','wt2654','wt0002','wt0001','wt0003'};
expr={'wt0001'};
% expr={'wt2654'};

% ----- wt0001 -----------------------------------
if ~isempty(strmatch('wt0001',expr))
  ci=ci+1;
  ri=0; 
  cd([WP.rootPath '\beta3_wtko\wt0001_04708']);
  for fi=1:4
    ri=ri+1; 
    AP=[]; DS=[];
    eval(['a0' sprintf('%.2i',fi) '_r1;']);
    AP_beta3_wtko;
    rmouse_APcheck;
    if isempty(ANPAR), ANPAR=AP; DSET=DS;
    else, ANPAR(ri,ci)=AP;  DSET(ri,ci)=DS;
    end
  end
  combine_tc(33.5);
end  

% ----- wt0002 -----------------------------------
if ~isempty(strmatch('wt0002',expr))
  ci=ci+1;
  ri=0; 
  cd([WP.rootPath '\beta3_wtko\wt0002_04707']);
  
  for fi=0:4
    ri=ri+1; 
    AP=[]; DS=[];
    eval(['a0' sprintf('%.2i',fi) '_r1;']);
    AP_beta3_wtko;
    rmouse_APcheck;
    if isempty(ANPAR), ANPAR=AP; DSET=DS;
    else, ANPAR(ri,ci)=AP;  DSET(ri,ci)=DS;
    end
  end
  combine_tc(63.6);
end

% ----- wt0003 -----------------------------------
if ~isempty(strmatch('wt0003',expr))
  ci=ci+1;
  ri=0; 
  cd([WP.rootPath '\beta3_wtko\wt0003_04730']);
  
  for fi=2:6
    ri=ri+1; 
    AP=[]; DS=[];
    eval(['a0' sprintf('%.2i',fi) '_r1;']);
    AP_beta3_wtko;
    rmouse_APcheck;
    if isempty(ANPAR), ANPAR=AP; DSET=DS;
    else, ANPAR(ri,ci)=AP;  DSET(ri,ci)=DS;
    end
  end
  combine_tc(65.3501);
end

% ----- wt2642 -----------------------------------
if ~isempty(strmatch('wt2642',expr))
  ci=ci+1;
  ri=0; 
  cd([WP.rootPath '\beta3_wtko\wt2642_03304']);
  
%  for fi=6:13    % all
  for fi=6:8      % all with behav scoring
    ri=ri+1; 
    AP=[]; DS=[];
    eval(['a0' sprintf('%.2i',fi) '_r1;']);
    AP_beta3_wtko;
    rmouse_APcheck;
    if isempty(ANPAR), ANPAR=AP; DSET=DS;
    else, ANPAR(ri,ci)=AP;  DSET(ri,ci)=DS;
    end
  end
  combine_tc(29+5/6);
end

% ----- wt2654 -----------------------------------
if ~isempty(strmatch('wt2654',expr))
  ci=ci+1;
  ri=0; 
  cd([WP.rootPath '\beta3_wtko\wt2654_03304']);
  
%  for fi=0:5    % all
  for fi=0:2     % all with behav scoring
%  for fi=0:1     % all with new behav scoring    
    ri=ri+1; 
    AP=[]; DS=[];
    eval(['a0' sprintf('%.2i',fi) '_r1;']);
    AP_beta3_wtko;
    rmouse_APcheck;
    if isempty(ANPAR), ANPAR=AP; DSET=DS;
    else, ANPAR(ri,ci)=AP;  DSET(ri,ci)=DS;
    end
  end
  combine_tc(36.33);
end

% ----- wt2008 -----------------------------------
% wt2008_02507 has no behavioral scoring
if ~isempty(strmatch('wt2008',expr))
  ci=ci+1;
  ri=0; 
  cd([WP.rootPath '\beta3_wtko\wt2008_02507']);
  for fi=7:10
    ri=ri+1; 
    AP=[]; DS=[];
    eval(['a0' sprintf('%.2i',fi) '_r1;']);
    AP_beta3_wtko;
    rmouse_APcheck;
    if isempty(ANPAR), ANPAR=AP; DSET=DS;
    else, ANPAR(ri,ci)=AP;  DSET(ri,ci)=DS;
    end
  end
  combine_tc(30.00); % 30 min is a guess!
end  

cd ..


% this section runs the combine_bv routine, generating a time course plot of behaviors
if size(ANPAR,2)>1
  combine_bv('mt',[33.5   63.6   65.3501   29+5/6   36.33]);
end