global DS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% animal's name/tag, experimental day and brief description of type of experiment - 
% should contain group-specific, unambiguous specifiers (like 'ko' or 'wt atropine') 
DS.aName=                      '2003-2-11; wt 2642';
% concentration of drug, if applied. special entries:
% NaN = no drug experiment
% 0   = control
% -1  = recovery
DS.conc=                       nan;
% path to files
DS.dpath=                      '\beta3_wtko\wt2642_03211';
% neuronal signals: abf file name without extension
DS.abfFn=                      '03211000';    
% some recording hardware (notably the Neuralynx) inverts signals by default.
% set this parameter to a nonzero value if that is the case
DS.rawSignalInverted=          1;
% Brief description(s) of experiment performed
DS.comment=                    'data are partly truncated in amplitude (max input range was +- 4 mV) ';
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
'Lynx1' ,0,  0,' ';...
'Lynx2' ,0,  100,'pyr';...   
'Lynx3' ,0,  200,' ';...    
'Lynx4' ,0,  300,' ';...   
'Lynx5' ,0, 400,' ';...  
'Lynx6' ,0, 500,'lm';...     
'Lynx7' ,0, 600,' ';...     
'Lynx8' ,0, 700,' ';...     
'Lynx9' ,0, 800,' ';...     
'Lynx10',0, 900,' ';...     
'Lynx11',0, 1000,' ';...     
'Lynx12',0, 1100,' ';...     
'Lynx13',0, 1200,' ';...     
'Lynx14',0, 1300,' ';...     
'Lynx15',0, 1400,' ';...     
'Lynx16',1, -inf,' '...
};

% The expected amplitude range of neuronal signals in mV. This value will 
% be used for scaling purposes (for the theta, gamma etc. streams). 
DS.nsRng=4;

% -------------------------- spare code -----------------------------------


