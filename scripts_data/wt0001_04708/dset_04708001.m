global AP DS
% animal's name/tag, experimental day and brief description of type of experiment - 
% should contain group-specific, unambiguous specifiers (like 'ko' or 'wt atropine') 
DS.aName=                      'wt #1; 2004-7-8; control';
% concentration of drug, if applied. special entries:
% NaN = no drug experiment
% 0   = control
% -1  = recovery
DS.conc=                       0;
% path to data files - either full path or relative path (in the latter case 
% rmouse_ini.m must contain the proper root path; the concatenation of both must be
% equivalent to the full path)
DS.dpath=                      '\beta3_wtko\wt0001_04708';
% neuronal signals: abf file name without extension
DS.abfFn=                      '04708001';    
% some recording hardware (notably the Neuralynx) inverts signals by default.
% set this parameter to a nonzero value if that is the case
DS.rawSignalInverted=          1;
% Comments
DS.comment=                    ' ';
% information about the NEURAL raw data channels that are OK in terms of data integrity: 
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
% <other>   - undefined or EEG 
% NOTE 1: the order of channels MUST be from dorsal to ventral
% NOTE 2: if there is a channel you are sure you'll never want to include in any kind of
% analysis, do not list it here
% NOTE 3: the distance between any two strata varies from animal to animal. So does, 
% as a consequence, the number of electrodes between those strata. Part of 
% the code combining results across animals will account for this variability by 
% linearly rescaling the interelectrode distances such that e.g. electrodes
% in s.pyr. are always at -1 and s.l.-m. electrodes always at 0 (arbitrary unit, 
% think of it as a 'functional depth' unit). Therefore, electrode location should be 
% determined with the greatest accuracy, and the strings must correspond EXACTLY to 
% those given above (e.g. sth. like 'pyr?' will not be recognized as the pyramidal layer). 
DS.rawCh={...
'IN 0' ,0,  0,' ';...
'IN 1' ,0,  100,' ';...   
'IN 2' ,0,  200,' ';...    
'IN 3' ,0,  300,' ';...   
'IN 4' ,0, 400,' ';...  
'IN 5' ,0, 500,' ';...     
'IN 6' ,0, 600,' ';...     
'IN 7' ,0, 700,' ';...     
'IN 8' ,0, 800,' ';...     
'IN 9' ,0, 900,' ';...     
'IN 10',0, 1000,' ';...     
'IN 11',0, 1100,'lm';...     
'IN 12',0, 1200,' ';...     
'IN 13',0, 1300,' ';...     
'IN 14',0, 1400,' ';...     
'IN 15',2, -inf,' '...
};

% The expected amplitude range of neuronal signals in mV. This value will 
% be used for scaling purposes (for the theta, gamma etc. streams). 
% ** NOTE: if one of the channels contains behavioral scoring trigger pulses 
%    DS.nsRng must cover the largest pulse amplitude!
DS.nsRng=10;
