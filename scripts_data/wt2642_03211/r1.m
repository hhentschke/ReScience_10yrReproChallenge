global AP

% invoke the data set-generating routine
dset_03211000;

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%  LENGTH OF DATA; CHANNELS; FILE NAMES; PRINT ETC
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% time block(s) to be processed in s (from t=0=start of recording).
% May be given as 'full'
AP.rawExcerpt=              'full';[0 200];

% the channels to be analyzed (must be a subset of DS.rawChan)
% NOTE: specify full names as given in the abf file header in cell 
% array of strings, like {'IN 0','IN 10'}
% AP.rawChAnNm=               DS.rawCh([1:10 12:16],1); % old
AP.rawChAnNm=               DS.rawCh([1:10 12],1);

% the 'principal' channel (=the one to which all others will be normalized) 
% Same rules as above apply
AP.rawChPrincNm=             {'Lynx6'};

% full path to results file(s) 
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
AP.bScoreFn=                '03211000_bscore';

% name of the behavioral scoring channel in abf file (will be ignored if text file 
% is used)
AP.rawChScoreNm=            {'niente'};

% s, point in time where behavioral recording was started, relative to the beginning
% of the neural recording (as indicated by timer in camera view). positive values
% mean the behavioral recording started after the neural recording.
AP.bbor=                     0;

% ms, the average (& assumed) delay between the onset of a behavior and the experimenter's 
% reaction (pressing a button)
AP.bsRt=                     300;

% types of scored behavior:
% column 1: type
% column 2: range of impulse voltages coding for this behavior
% column 3: color (used in plots), either single char or rgb arr
% column 4: symbol type to be used in plots
% NOTES: 
% - 'bad' behavior = segments you wish to be omitted from analysis for 
%   whatever reasons
% - impulse voltages in column 2: currently, the upper limit is ignored for detection 
%   of impulses in abf files (but not in txt files); make sure the lower limit is
%   correct
AP.behavType={...
'grooming',  [1.0 2.4],  [.5 .4 1],   '^';...
'exploring', [6     9],  [.5 1 .5],   'o';...
'immobile',  [2.5   4],  [.6 .4 .1],  's';...
'bad',       [-1   -9],  [.9 0 0],    '+'...
};

% Hz, corner frequency of highpass filter to apply to behavioral scoring data 
% (in some instances the baseline drifts slowly). Set to [] to skip filtering.
AP.bsCFreq=                  .1;

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%        DETECTION OF ARTIFACTS
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% the neuronal channels to be monitored for periods of bad signal quality (noise) 
% and/or other unwanted periods. Must be a subset of AP.rawChAnNm
AP.rawChMonNm=               {'Lynx6','Lynx12'};

% mV, thresholds for detection of artifacts on raw channels specified in 
% AP.rawChMonNm 
AP.afThresh= -1*                  [3 -3.8];

% Hz, corner frequencies of lowpass filter to apply to channel monitored for artifacts
% set to [] if no filter shall be applied
AP.afCFreq=                   [];

% seconds, the interval around an artifactual event/period (noise) to omit from analysis
AP.afWin=                      [-.25 .25];

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% FREQUENCIES & FILTERS
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ------ center frequency for butterworth bandstop filter. If nonempty, the gamma stream 
% will be filtered at corner frequencies +/-1 this value with high rolloff (set to [] 
% to prevent bandstop filtering)
AP.lineFreq=                   [];

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% SPECTRAL & CORRELATION ANALYSIS PARAMETERS
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
