% collect all dsets and ANPARs 
global ANPAR DSET AP DS WP


% ANPAR and DSET will be set up as 2D or 3D struct arrays. All data sets to
% be averaged will reside in rows. Different rows will code for
% concentrations, in this order: 
%                   control (no drug)
%                   drug
%                   recovery
%
% Different 'slices' account for different genotypes in this order:
%                   wild type
%                   knockout/in


ANPAR=[];
DSET=[];

% learn data root path from rmouse ini routine
rmouse_ini;
i3=0;

% --------------- WT --------------------------------
% --------------- WT --------------------------------
ri=1; 
ci=0;
i3=i3+1;

ci=ci+1;
cd([WP.rootPath '\beta3_wtko\wt0001_04708']);
AP=[]; DS=[];
AP_beta3_wtko;
% a001_r1;
a001_exc1;
rmouse_apcheck;
if isempty(ANPAR), ANPAR=AP; DSET=DS;
else, ANPAR(ri,ci,i3)=AP;  DSET(ri,ci,i3)=DS;
end

ci=ci+1;
cd([WP.rootPath '\beta3_wtko\wt0002_04707']);
AP=[]; DS=[];
AP_beta3_wtko;
a000_r1;
rmouse_apcheck;
if isempty(ANPAR), ANPAR=AP; DSET=DS;
else, ANPAR(ri,ci,i3)=AP;  DSET(ri,ci,i3)=DS;
end

ci=ci+1;
cd([WP.rootPath '\beta3_wtko\wt0003_04730']);
AP=[]; DS=[];
AP_beta3_wtko;
a003_r1;
rmouse_apcheck;
if isempty(ANPAR), ANPAR=AP; DSET=DS;
else, ANPAR(ri,ci,i3)=AP;  DSET(ri,ci,i3)=DS;
end

ci=ci+1;
cd([WP.rootPath '\beta3_wtko\wt2642_03211']);
AP=[]; DS=[];
AP_beta3_wtko;
r1;
rmouse_apcheck;
if isempty(ANPAR), ANPAR=AP; DSET=DS;
else, ANPAR(ri,ci,i3)=AP;  DSET(ri,ci,i3)=DS;
end

ci=ci+1;
cd([WP.rootPath '\beta3_wtko\wt2654_03305']);
AP=[]; DS=[];
AP_beta3_wtko;
r1;
rmouse_apcheck;
if isempty(ANPAR), ANPAR=AP; DSET=DS;
else, ANPAR(ri,ci,i3)=AP;  DSET(ri,ci,i3)=DS;
end

ci=ci+1;
cd([WP.rootPath '\beta3_wtko\wt2008_02220']);
AP=[]; DS=[];
AP_beta3_wtko;
r1;
rmouse_apcheck;
if isempty(ANPAR), ANPAR=AP; DSET=DS;
else, ANPAR(ri,ci,i3)=AP;  DSET(ri,ci,i3)=DS;
end

ci=ci+1;
cd([WP.rootPath '\beta3_wtko\wt2009_02319']);
AP=[]; DS=[];
AP_beta3_wtko;
r1;
rmouse_apcheck;
if isempty(ANPAR), ANPAR=AP; DSET=DS;
else, ANPAR(ri,ci,i3)=AP;  DSET(ri,ci,i3)=DS;
end

% the bad ones: 
% cd d:\rmouse\wt2141_01d18
% cd d:\rmouse\wt2141_02510
% cd d:\rmouse\wt2001_02424
% cd d:\rmouse\wt2206_02725

% --------------- KO --------------------------------
% --------------- KO --------------------------------
ri=1; 
ci=0;
i3=i3+1;

ci=ci+1;
cd([WP.rootPath '\beta3_wtko\ko2472_02o11']);
AP=[]; DS=[];
AP_beta3_wtko;
a000_r1;
rmouse_apcheck;
if isempty(ANPAR), ANPAR=AP; DSET=DS;
else, ANPAR(ri,ci,i3)=AP;  DSET(ri,ci,i3)=DS;
end

ci=ci+1;
cd([WP.rootPath '\beta3_wtko\ko2484_02927']);
AP=[]; DS=[];
AP_beta3_wtko;
a002_r1;
rmouse_apcheck;
if isempty(ANPAR), ANPAR=AP; DSET=DS;
else, ANPAR(ri,ci,i3)=AP;  DSET(ri,ci,i3)=DS;
end

ci=ci+1;
cd([WP.rootPath '\beta3_wtko\ko3382_02801']);
AP=[]; DS=[];
AP_beta3_wtko;
r1;
rmouse_apcheck;
if isempty(ANPAR), ANPAR=AP; DSET=DS;
else, ANPAR(ri,ci,i3)=AP;  DSET(ri,ci,i3)=DS;
end

ci=ci+1;
cd([WP.rootPath '\beta3_wtko\ko2483_03303']);
AP=[]; DS=[];
AP_beta3_wtko;
r1;
rmouse_apcheck;
if isempty(ANPAR), ANPAR=AP; DSET=DS;
else, ANPAR(ri,ci,i3)=AP;  DSET(ri,ci,i3)=DS;
end

ci=ci+1;
cd([WP.rootPath '\beta3_wtko\ko3363_02313']);
AP=[]; DS=[];
AP_beta3_wtko;
r1;
rmouse_apcheck;
if isempty(ANPAR), ANPAR=AP; DSET=DS;
else, ANPAR(ri,ci,i3)=AP;  DSET(ri,ci,i3)=DS;
end

ci=ci+1;
cd([WP.rootPath '\beta3_wtko\ko2445_03214']);
AP=[]; DS=[];
AP_beta3_wtko;
r1;
rmouse_apcheck;
if isempty(ANPAR), ANPAR=AP; DSET=DS;
else, ANPAR(ri,ci,i3)=AP;  DSET(ri,ci,i3)=DS;
end

ci=ci+1;
cd([WP.rootPath '\beta3_wtko\ko3362_02528']);
AP=[]; DS=[];
AP_beta3_wtko;
r1;
rmouse_apcheck;
if isempty(ANPAR), ANPAR=AP; DSET=DS;
else, ANPAR(ri,ci,i3)=AP;  DSET(ri,ci,i3)=DS;
end


cd ..



% -------------- plausibility checks --------------------------------------
% check whether first conc is 0 
conc=cat(2,DSET(1,:).conc);
if any(conc), error('check control concentrations (DSET.conc)'); end 

% if any two combinations of abf file name AND concentration are identical
% we are very likely to have a foul data set
ct=1;
for g=1:numel(DSET)
  if ~isempty(DSET(g).abfFn)
    eNm{ct}=[DSET(g).abfFn ' ' sprintf('%3.1f',DSET(g).conc)]; 
    ct=ct+1;
  end
end

if length(eNm)~=length(unique(eNm))
  warndlg('There are at least two identical experiments in the current collection of data (as judged by DS.aName)');
  strvcat(eNm)
end

% ------------------- RInfo ---------------------------------
% NOTE:
% - RInfo is a required input into combine_r (which collects data from all
%   experiments listed in DSET/ANPAR and saves them as a cell array in a 
%   file) 
% - the factors below ('RInfo.name') differ slightly in how they influence
%   collection of data and/or statistics
% - 'behavior': any behavioral level not listed below (RInfo(2).level) will
%   not be collected by combine_r (because there are only two out of seven
%   or more which will ever be interesting, so there is no point in 
%   collecting the others). 
% - 'drug': the structure of DSET/ANPAR reflects the number of levels 
%   (=drug treatments); inconsistencies will cause an error
% 'genotype': same as 'drug'

% behavior - choose the levels you want to distil
RInfo(1).name='behavior';
RInfo(1).level={'immobile';'exploring'};

RInfo(2).name='drug'; 
RInfo(2).level={'control'};

RInfo(3).name='genotype';
RInfo(3).level={'WT','KO'};

% § check what happens if number of wt and ki diverges
RInfo(4).name='subject'; 
RInfo(4).level={DSET(1,:,:).aName};


% --- concluding code checking for congruence between RInfo and DSET/ANPAR --
tmpIx=strmatch('drug',{RInfo.name});
if ~isempty(tmpIx)
  if length(RInfo(tmpIx).level)~=size(DSET,1),
    warndlg('the number of rows of DSET/ANPAR does not match the number of drug levels specified in RInfo');
  end
end
% § extend to 3D
tmpIx=strmatch('subject',{RInfo.name});
if ~isempty(tmpIx)
  if numel(RInfo(tmpIx).level)~=numel(DSET(1,:,:)),
    warndlg('the number of columns of DSET/ANPAR does not match the number of animals specified in RInfo');
  end
end
tmpIx=strmatch('genotype',{RInfo.name});
if ~isempty(tmpIx)
  if length(RInfo(tmpIx).level)~=size(DSET,3),
    warndlg('the number of ''slices'' of DSET/ANPAR does not match the number of genotypes specified in RInfo');
  end
end

% -------------------- cleanup! -------------------------------------------
clear pool ri ci i3 expr* nWt nKi nAll nIncr conc eNm apf* cur* g subDi cOrd proceed fInd ct

