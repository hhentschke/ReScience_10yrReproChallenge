global AP

% invoke the data set-generating routine
dset_04708001;

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%  LENGTH OF DATA; CHANNELS; FILE NAMES; PRINT ETC
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% time block(s) to be processed in s (from t=0=start of recording).
% May be given as 'full'
AP.rawExcerpt=              [2.5 3.5]*60; % mixed behavior
% AP.rawExcerpt=              [0 15]; % mixed behavior
% AP.rawExcerpt=              [0 2]*60; % only exploring
% AP.rawExcerpt=              [171.3 172.8]; %


% the channels to be analyzed (must be a subset of DS.rawChan)
% NOTE: specify full names as given in the abf file header in cell 
% array of strings, like {'IN 0','IN 10'}
% ** IN 3 contains too much 60 Hz pickup, IN 15 is scoring channel
AP.rawChAnNm=               DS.rawCh([1:3 5:12],1);

% the 'principal' channel (serving as a reference for e.g. CC computations) 
% Same rules as above apply
AP.rawChPrincNm=             {'IN 11'};


% path to results file(s) - either full path or relative path (in the latter case 
% rmouse_ini.m must contain the proper root path; the concatenation of both must be
% equivalent to the full path)
AP.resPath=                  DS.dpath;

% name of the results file(s) 
AP.resFn=                    [DS.abfFn '_proc_test'];

% directory for intermediate results files - will be generated automatically 
% if not existent. make sure there is enough disk space (at least ten times 
% the size of one experimental file)
AP.strmDir=                  [DS.dpath '\streams'];

% name of the log file (which collects information about each run of rmouse.m
% on the data)
AP.logFn=                    [DS.abfFn '_logtest.txt'];


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
AP.bScoreFn=                '2005-02-28_0002';

% name of the behavioral scoring channel in abf file (will be ignored if text file 
% is used)
AP.rawChScoreNm=            {'score'};

% s, point in time where behavioral recording was started, relative to the beginning
% of the neural recording (as indicated by timer in camera view). positive values
% mean the behavioral recording started after the neural recording.
AP.bbor=                     0;

% ms, the assumed average delay between the onset of a behavior and the experimenter's 
% reaction (pressing a button)
AP.bsRt=                       400;

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

% Hz, corner frequency of highpass filter to apply to behavioral scoring data 
% (in some instances the baseline drifts slowly). Set to [] to skip filtering.
AP.bsCFreq=                  [];

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%        DETECTION OF ARTIFACTS
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% the neuronal channel(s) to be monitored for periods of bad signal quality (noise) 
% and/or other unwanted periods. Must be a subset of AP.rawChAnNm
AP.rawChMonNm=               {'IN 11'};

% mV, threshold(s) for detection of artifacts on raw channels specified in 
% AP.rawChMonNm 
AP.afThresh= -1*                  [-2];

% Hz, corner frequencies of lowpass filter to apply to channel monitored for artifacts
% set to [] if no filter shall be applied
AP.afCFreq=                   [];

% seconds, the interval around an artifactual event/period (noise) to omit from analysis
AP.afWin=                      [-1 1];

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% FREQUENCIES & FILTERS
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ------ center frequency for butterworth bandstop filter. If nonempty, the gamma stream 
% will be filtered at corner frequencies +/-1 this value with high rolloff (set to [] 
% to prevent bandstop filtering)
AP.lineFreq=                   [60];

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% SPECTRAL & CORRELATION ANALYSIS PARAMETERS
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



% this is a partial AP (Analysis Parameter) set to be used in batch processing 

AP.job={...
  'det_artifacts',...       % detect artifacts in neural signals & generate time stamp list thereof
  'gen_btsl',...            % generate behavioral time stamp list
  'chechao',...             % check proper order of channels by crosscorrelating raw data
  'filter&hilbert',...      % generate the various streams (theta, gamma, etc.)
  'gen_avDiff',...          % generate a stream (essentially, slopes averaged across channels) useful for detetction of epileptiform activity
  'thetaPeaks',...          % compute amplitude & time of occurrence of peaks in theta stream
  'gammaEnvPeaks',...       % compute amplitude & time of occurrence of peaks in gammaEnv stream  
  'specgram',...            % compute & plot spectrogram (time course of spectral power)
  'seg_thetaPeakReg',...    % compute theta regularity based on variance of peak amplitude and inter-peak interval
  'seg_gammaEnvPeakReg',... % compute gammaEnv regularity based on variance of peak amplitude and inter-peak interval
  'seg_rawSpecPow',...      % compute spectra & derived parameters (all combinations of channels)
  'seg_gammaEnvSpecPow',... % compute spectra & derived parameters for gamma envelope (all combinations of channels)  
  'seg_gammaNarrowEnvSpecPow',... % compute spectra & derived parameters for narrow gamma envelope (all combinations of channels)  
  'seg_deltaCC',...         % compute delta crosscorrelations (all combinations of channels)
  'seg_thetaCC',...         % compute theta crosscorrelations (all combinations of channels)
  'seg_thetaLoEnvCC',...    % compute theta_lo envelope crosscorrelations (all combinations of channels)
  'seg_thetaHiEnvCC',...    % compute theta_hi envelope crosscorrelations (all combinations of channels)
  'seg_gammaCC',...         % compute gamma crosscorrelations (all combinations of channels)
  'seg_gammaEnvCC',...      % compute gamma envelope crosscorrelations (all combinations of channels)
  'seg_gammaNarrowCC',...   % compute gamma narrow crosscorrelations (all combinations of channels)
  'seg_gammaNarrowEnvCC',...% compute gamma narrow envelope crosscorrelations (all combinations of channels)
  'seg_thetaGammaEnvCC',... % compute crosscorrelations between theta and gamma envelope (for each channel)
  'seg_deltaThetaEnvCC',... % compute crosscorrelations between delta and theta envelope (for each channel)
  'seg_rawGammaEnvCoh',...  % compute coherence between raw signal and gamma envelope (all combinations of channels)
  'seg_cc_ms2rad',...       % convert theta-related phase lags from ms to radians, creating new results fields (NOTE: will be called automatically if any of 'seg_rawSpecPow','seg_thetaCC','seg_gammaEnvCC','seg_thetaGammaEnvCC' is called)
  'sumFig',...              % produce summary plots
  'rawPlot',...             % plot of raw data excerpts
  'ovRawPlot',...           % overview plot of raw data (the whole excerpt to be analyzed)
  'csdPlot',...             % contour plot of current source densities from theta stream excerpts  
  'tcFig',...               % plots time course of segment-wise computed, selected parameters 
  'tcPrincComp',...         % computes & plots principal components from tupels of segment-wise computed, selected parameters 
  'tcThetaPeak',...         % produce blob-plots depicting time course of theta peaks (the plots have an aesthetic appeal and may reveal variability of theta)
};

AP.job={...
  'gen_btsl',...            % generate behavioral time stamp list
  'seg_thetaCC',...         % compute theta crosscorrelations (all combinations of channels)
  };

% determines how to deal with previously computed results of the 'seg_' jobs:
% 'replace' means they will be discarded entirely and replaced by the results 
%   of the current run. This should be the default.
% 'append' means that the results of the current run will be appended to or 
%   overwrite results of previous runs. 
% *** see NOTE of AP.job above! ***
AP.saveMode=                'append';
AP.saveMode=                'replace';

% save plot(s) to file in what format? '-dpsc2' (color postscript level 2) is
% strongly recommended. multiple formats can be specified; in this case they must 
% be listed in a cell array, like {'-dpsc2','-djpeg90'}. set to [] if plots 
% shall not be saved.
AP.printas=                  {'-djpeg90'};
AP.printas=                  [];

% the significance of the CC between theta and gamma envelope is evaluated by
% computing the CC between theta and phase-shuffled versions of the gamma envelope.
% The parameter below specifies on how many shuffled versions of gammaEnv this
% evaluation should be based. Set to zero to skip shuffling.
AP.ccNShuffle=                  0;   
