
% ------------------- statistics settings ---------------------------------
% NOTE:
% - idv is a required input into rdeal (which does statistics and plots)
% - the factors below ('idv.name') differ slightly in how they influence
%   collection of data and/or statistics
% - 'rec site': no matter which levels are listed, combine_r will gather 
%   ALL that there are in the data; it is in rdeal where the collected data
%   will be pruned. This is so because statistical analysis strategies
%   depend strongly on the set of available rec sites, which in turn
%   depends on the pool of animals investigated, so we need maximal
%   flexibility here. The number of available levels of both 'rec site' and 
%   'behavior' cannot be inferred from the structure of DSET/ANPAR because
%   they are factors within data sets. 
% - 'behavior': any behavioral level not listed below (idv(2).level) will
%   not be collected by combine_r (because there are only two out of seven
%   or more which will ever be interesting, so there is no point in 
%   collecting the others). 
% - 'drug': similar to 'rec site', combine_r will collect all and rdeal
%   will prune, which is more flexible because collection of data with 
%   combine_r takes awfully long while pruning is faster than an eye blink.
%   However, in contrast to 'rec site', the structure of DSET/ANPAR
%   reflects the number of levels (=drug treatments); inconsistencies will
%   cause an error
% 'genotype': same as 'drug'

idv(1).name='rec site';
% rec sites; conventions: depth in mm; 0 is principal channel, spacing 0.1,
% values decreasing in dorsal direction
idv(1).level={-0.6; -0.5; -0.4; -0.3; -0.2; -0.1; 0.0};
% idv(1).level={-0.5; -0.4; -0.3; -0.2; -0.1; 0.0};
idv(1).pCol={[]};
idv(1).pSymb={[]};
idv(1).treat=repmat({'stats'},size(idv(1).level));

% 2nd factor behavior - choose the levels you want to distil
idv(2).name='behavior';
idv(2).level={'immobile';'exploring'};
idv(2).pCol={[.6 .4 .1];[.5 1 .5]};
idv(2).pSymb={'s';'o'};
idv(2).treat={'stats';'stats'};
idv(2).treat={'separate';'separate'};



% third factor
idv(3).name='drug'; 
idv(3).level={'control'};
idv(3).pCol={'k'};
idv(3).pSymb={'o'};
idv(3).treat={'ignore'};

% fourth factor
idv(4).name='genotype';
idv(4).level={'WT','KO'};
idv(4).pCol={'k',[.6 .6 .6]};
idv(4).pSymb={'o','s'};
idv(4).treat={'stats';'stats'};

% the fifth factor, the individuals 
idv(5).name='subject'; 
idv(5).level={DSET(1,:,:).aName};
idv(5).pCol={[]};
idv(5).pSymb={[]};
idv(5).treat=repmat({'ignore'},size(idv(5).level));

% --- concluding code checking for congruence between idv and DSET/ANPAR --
tmpIx=strmatch('drug',{idv.name});
if ~isempty(tmpIx)
  if length(idv(tmpIx).level)~=size(DSET,1),
    warndlg('the number of rows of DSET/ANPAR does not match the number of drug levels specified in idv');
  end
end

% -------------------- cleanup! -------------------------------------------
clear pool ri ci i3 expr* nWt nKi nAll nIncr conc eNm apf* cur* g subDi cOrd proceed fInd ct
