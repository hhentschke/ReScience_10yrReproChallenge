global AP

% invoke the data set-generating routine
dset_04707002;


% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%  LENGTH OF DATA; CHANNELS; FILE NAMES; PRINT ETC
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% time block(s) to be processed in s (from t=0=start of recording).
% May be given as 'full'
AP.rawExcerpt=              'full';[0 400];

% the channels to be analyzed (must be a subset of DS.rawChan)
% NOTE: specify full names as given in the abf file header in cell 
% array of strings, like {'IN 0','IN 10'}
% ** IN 3 contains too much 60 Hz pickup, IN 15 is scoring channel
AP.rawChAnNm=               DS.rawCh([1:3 5:end-1],1);

% the 'principal' channel (serving as a reference for e.g. CC computations) 
% Same rules as above apply
AP.rawChPrincNm=             {'IN 10'};


% path to results file(s) - either full path or relative path (in the latter case 
% rmouse_ini.m must contain the proper root path; the concatenation of both must be
% equivalent to the full path)
AP.resPath=                  DS.dpath;

% name of the results file(s) 
AP.resFn=                    [DS.abfFn '_proc_r1'];

% directory for intermediate results files - will be generated automatically 
% if not existent. make sure there is enough disk space (at least ten times 
% the size of one experimental file)
AP.strmDir=                  [DS.dpath '\streams'];

% name of the log file (which collects information about each run of rmouse.m
% on the data)
AP.logFn=                    [DS.abfFn '_log.txt'];

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%        BEHAVIOR
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% data for behavioral scoring: file name without extension 
% two data formats are supported:
% 1. abf, in which case the file must contain trigger pulses as described below
% 2. txt, a text file which must contain in its first column time stamps (in ms) 
%    and in its second column a number corresponding to the trigger level (listed
%    in column 2 of AP.bScoreFn below). The time stamps need not be sorted.
% NOTES
% a) If two files with identical names but different extensions are found the 
%    txt file is used
% b) In cases where behavioral scoring data was recorded in parallel with and 
%    in the same file as neural data, that file name must be specified here again
AP.bScoreFn=                '04707002';

% name of the behavioral scoring channel in abf file (will be ignored if text file 
% is used)
AP.rawChScoreNm=            {'IN 15'};

% s, point in time where behavioral recording was started, relative to the beginning
% of the neural recording (as indicated by timer in camera view). positive values
% mean the behavioral recording started after the neural recording.
AP.bbor=                     0;

% ms, the assumed average delay between the onset of a behavior and the experimenter's 
% reaction (pressing a button)
AP.bsRt=                       250;

% types of scored behavior:
% column 1: type
% column 2: range of impulse voltages coding for this behavior
% column 3: color (used in plots), either single char or rgb array
% column 4: symbol to be used in plots
% NOTES: 
% - 'bad' behavior = segments you wish to be omitted from analysis for 
%   whatever reasons
% - impulse voltages in column 2: currently, the upper limit is ignored for detection 
%   of impulses in abf files (but not in txt files); make sure the lower limit is
%   correct
AP.behavType={...
    'grooming',  [1.0 3],  [.5 .4 1],  '^';...
    'exploring', [7  11],  [.5 1 .5],   'o';...
    'immobile',  [3.5 6],  [.6 .4 .1],  's';...
    'bad',       [-1 -11],  [.9 0 0],    '+'...
};

% ms, the trigger signal's rise time from base line to threshold (in behavioral 
% scoring data) ==>> make 100% sure this is the slowest rise time the trigger pulses 
% will ever show, otherwise categorization of segments will get messed up seriously
AP.trigTau=                  20;

% Hz, corner frequency of highpass filter to apply to behavioral scoring data 
% (in some instances the baseline drifts slowly). Set to [] to skip filtering.
AP.bsCFreq=                  [];

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%        DETECTION OF ARTIFACTS
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% the neuronal channel(s) to be monitored for periods of bad signal quality (noise) 
% and/or other unwanted periods. Must be a subset of AP.rawChAnNm
AP.rawChMonNm=               {'IN 14'};

% mV, threshold(s) for detection of artifacts on raw channels specified in 
% AP.rawChMonNm 
AP.afThresh=                  [-2];

% Hz, corner frequencies of lowpass filter to apply to channel monitored for artifacts
% set to [] if no filter shall be applied
AP.afCFreq=                   [];

% seconds, the interval around an artifactual event/period (noise) to omit from analysis
AP.afWin=                      [-1 1];

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% FREQUENCIES & FILTERS
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% SPECTRAL & CORRELATION ANALYSIS PARAMETERS
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
