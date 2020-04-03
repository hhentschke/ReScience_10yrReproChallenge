% ************************************************************************************
%                  This dset is the master template for other dsets                 
% ************************************************************************************
% - it contains ALL fields of DS
% - it contains the most comprehensive and up-to-date information (including
%   examples of values) about the fields
% - all fields are set to {[]} (cell arrays) or [] (other data types) in the bottom
%   section 
% - it is recommended to run it before invoking any DS, which ensures that all 
%   fields of DS at least exist (which avoids errors when DS are concatenated).  
% ************************************************************************************

global DS AP

% set up template DS
if ~strcmpi(mfilename,'dset_template'), DSET_template; end

% animal's name/tag, experimental day and brief description of type of experiment - 
% should contain group-specific, unambiguous specifiers (like 'ko' or 'wt atropine') 
DS.aName=                      'wt 2001; 2002-4-24; control';
% concentration of drug, if applied. special entries:
% NaN = no drug experiment
% 0   = control
% -1  = recovery
DS.conc=                       nan;
% path to data files - either full path or relative path (in the latter case 
% rmouse_ini.m must contain the proper root path; the concatenation of both must be
% equivalent to the full path)
DS.dpath=                      '\beta3_wtko\wt2001_02424';
% neuronal signals: abf file name without extension
DS.abfFn=                      '02424000';    
% some recording hardware (notably the Neuralynx) inverts signals by default.
% set this parameter to a nonzero value if that is the case
DS.rawSignalInverted=          1;
% Comments
DS.comment=                    ' ';
% information about the NEURAL raw data channels:
% column 1: channel name as given in abf file header
% column 2: type of data (0=LFP, 1=EEG, 2=other (..to be refined))
% column 3: location (distance relative to arbitrary reference in um) 
% column 4: location: name of stratum. Must be one of
% 'cx'      - neocortex
% 'oriens'  - s. oriens or alveus
% 'pyr'     - s. pyramidale
% 'rad'     - s.radiatum
% 'lm'      - s. lacunosum-moleculare
% 'mol'     - DG, molecular layer
% 'gran'    - DG, granular layer
% 'polym'   - DG, polymorphic layer (hilus)
%  other    - undefined or EEG 
%        
%        ******* please read and make sure you understand notes 1 & 2 below!! ******
%
% NOTE 1: the order of channels MUST be from dorsal to ventral
% NOTE 2: even if there is a LFP channel you are sure you'll never want to
% include in any kind of analysis because its signals are bad, list it
% here. You can exclude it from analysis in AP.rawChAnNm. The reason lies
% in the way rmouse organizes results data according to electrode depth
% NOTE 3: the distance between any two strata varies from animal to animal.
% So does, as a consequence, the number of electrodes between those strata.
% Future versions of rmouse could possibly account for this variability by 
% linearly rescaling the interelectrode distances such that e.g. electrodes
% in s.pyr. are always at -1 and s.l.-m. electrodes always at 0 (arbitrary
% unit, think of it as a 'functional depth' unit). Therefore, if
% information about electrode location is available, it may be listed in
% the fourth column. If you do so, electrode location should be determined
% with the greatest accuracy, and the strings must correspond EXACTLY to 
% those given above (e.g. sth. like 'pyr?' will not be recognized as the
% pyramidal layer). Currently, however, rmouse does not use this
% information.

% Here's one frequently used template 
DS.rawCh={...
'IN 0' ,0,  0,' ';...
'IN 1' ,0,  100,' ';...   
'IN 2' ,0,  200,' ';...    
'IN 3' ,0,  300,' ';...   
'IN 4' ,0, 400,' ';...  
'IN 5' ,0, 500,' ';...     
'IN 6' ,0, 600,' ';...     
'IN 7' ,0, 700,' ';...     
'IN 8' ,0, 800,'cx';...     
'IN 9' ,0, 900,'pyr';...     
'IN 10',0, 1000,' ';...     
'IN 11',0, 1100,'lm';...     
'IN 12',0, 1200,'lm?';...     
'IN 13',0, 1300,'mol?';...     
'IN 14',0, 1400,'mol?';...     
'IN 15',1, -inf,' '...
};

% ..and another
DS.rawCh={...
'Lynx1' ,0,  0,' ';...
'Lynx2' ,0,  100,' ';...   
'Lynx3' ,0,  200,' ';...    
'Lynx4' ,0, 300,' ';...   
'Lynx5' ,0, 400,'pyr';...  
'Lynx6' ,0, 500,'o';...     
'Lynx7' ,0, 600,'r';...     
'Lynx8' ,0, 700,'r';...     
'Lynx9' ,0, 800,'lm';...     
'Lynx10',0, 900,' ';...     
'Lynx11',0, 1000,' ';...     
'Lynx12',0, 1100,' ';...     
'Lynx13',0, 1200,' ';...     
'Lynx14',0, 1300,' ';...     
'Lynx15',0, 1400,' ';...     
'Lynx16',1, 1500,' '...
};

% The expected amplitude range of signals in mV. This value will be used 
% for scaling purposes (generation and upload of theta, gamma etc. streams). 
% ** NOTE: make sure that the range is wide enough to cover all neurological
% signals that shall be analyzed. Behavioral scoring trigger pulses 
% (which are usually much larger than neuronal signals) do not count
% because they will not be converted. Neither should the range be so wide
% as to encompass artifacts - these will be ignored in the analysis anyways.
DS.nsRng=4;

% ---------------------------------------------------------------------------
% ***** IN ANY DS TO BE USED FOR REAL THE CODE BELOW MAY BE ENTIRELY DELETED 
% OR MUST BE LEFT UNALTERED *****
if strcmpi(mfilename,'dset_template')
  s=fieldnames(DS);
  for ns=1:length(s)
    eval(['c=iscell(DS.' s{ns} ');']);
    if c, eval(['DS.' s{ns} '={[]};']);
    else eval(['DS.' s{ns} '=[];']);
    end
  end
  clear s ns c
end