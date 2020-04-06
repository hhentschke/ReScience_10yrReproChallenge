% collect all dsets and ANPARs 
global ANPAR DSET AP DS

% ANPAR and DSET will be set up as 2D or 3D struct arrays. All data sets to
% be averaged will reside in rows. Different rows will code for
% concentrations, in this order: 
%                   control (no drug)
%                   drug


ANPAR=[];
DSET=[];

% script root path 
script_root_path = "d:\hh\projects_programming\ten-years\ReScience_10yrReproChallenge\scripts_data";

ci=0;

expName = "wt0001_04708";
% reset row index, increase column counter for new subject
ri=0; 
ci = ci + 1;

% increase row counter for new drug condition
ri = ri + 1;
AP=[]; 
DS=[];
AP_beta3_wtko;
run(fullfile(script_root_path, expName, "a001_exc1.m"));
checkin_apds(ri, ci);

ri = ri + 1;
AP=[]; 
DS=[];
AP_beta3_wtko;
run(fullfile(script_root_path, expName, "a003_r1.m"));
checkin_apds(ri, ci);


expName = "wt0002_04707";
% reset row index, increase column counter for new subject
ri=0; 
ci = ci + 1;
% increase row counter for new drug condition
ri = ri + 1;
AP=[]; 
DS=[];
AP_beta3_wtko;
run(fullfile(script_root_path, expName, "a001_r1.m"));
checkin_apds(ri, ci);

ri = ri + 1;
AP=[]; 
DS=[];
AP_beta3_wtko;
run(fullfile(script_root_path, expName, "a002_r1.m"));
checkin_apds(ri, ci);


expName = "wt0003_04730";
% reset row index, increase column counter for new subject
ri=0; 
ci = ci + 1;
% increase row counter for new drug condition
ri = ri + 1;
AP=[]; 
DS=[];
AP_beta3_wtko;
run(fullfile(script_root_path, expName, "a003_r1.m"));
checkin_apds(ri, ci);

ri = ri + 1;
AP=[]; 
DS=[];
AP_beta3_wtko;
run(fullfile(script_root_path, expName, "a005_r1.m"));
checkin_apds(ri, ci);


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
RInfo(2).level={'control', 'atropine'};

RInfo(3).name='genotype';
RInfo(3).level={'WT'};

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

% cleanup 
clear ci conc ct eNm expName g ri script_root_path tmpIx
clear global DS AP

disp(mfilename + " finished successfully:")
whos 


function checkin_apds(ri, ci)
global ANPAR DSET AP DS
rmouse_apcheck;
if isempty(ANPAR)
    ANPAR=AP; 
    DSET=DS;
else
    ANPAR(ri,ci)=AP;  
    DSET(ri,ci)=DS;
end
end
