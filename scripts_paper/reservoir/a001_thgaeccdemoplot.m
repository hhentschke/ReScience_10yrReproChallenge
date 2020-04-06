global AP

% invoke the data set-generating routine
dset_04708001;

AP.job={...
  'det_artifacts',...       % detect artifacts in neural signals & generate time stamp list
  'gen_btsl',...            % generate behavioral time stamp list
};

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%  LENGTH OF DATA; CHANNELS; FILE NAMES; PRINT ETC
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% time block(s) to be processed in s (from t=0=start of recording).
% May be given as 'full'
AP.rawExcerpt=              'full';[0 30];

% the channels to be analyzed (must be a subset of DS.rawChan)
% NOTE: specify full names as given in the abf file header in cell 
% array of strings, like {'IN 0','IN 10'}
AP.rawChAnNm=               {'IN 11'};

% the 'principal' channel (serving as a reference for e.g. CC computations) 
% Same rules as above apply
AP.rawChPrincNm=             {'IN 11'};

% determines how to deal with previously computed results of the 'seg_' jobs:
% 'replace' means they will be discarded entirely and replaced by the results 
%   of the current run. This should be the default.
% 'append' means that the results of the current run will be appended to or 
%   overwrite results of previous runs. 
% *** see NOTE of AP.job above! ***
AP.saveMode=                 'replace';'append';

% path to results file(s) - either full path or relative path (in the latter case 
% rmouse_ini.m must contain the proper root path; the concatenation of both must be
% equivalent to the full path)
AP.resPath=                  'd:\projects\madison\rmouse\paper_atropine\rawFig';
AP.resPath=                  'd:\hh';

% name of the results file(s) 
AP.resFn=                    [DS.abfFn '_thgaeCCDemo'];

% directory for intermediate results files - will be generated automatically 
% if not existent. make sure there is enough disk space (at least ten times 
% the size of one experimental file)
AP.strmDir=                  [DS.dpath '\streams'];

% name of the log file (which collects information about each run of rmouse.m
% on the data)
AP.logFn=                    [DS.abfFn '_rubble.txt'];

% save plot(s) to file in what format? '-dpsc2' (color postscript level 2) is
% strongly recommended for high-quality figures, whereas '-djpeg90' (jpeg with 
% quality level of 90) is more efficient for everyday use. Multiple formats can 
% be specified; in this case they must be listed in a cell array, like 
% {'-dpsc2','-djpeg90'}. Set to [] if plots shall not be saved.
AP.printas=                  {'-dpsc2'};

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
% ------ definition of frequency bands (Hz)
AP.theta=                      [4 12];
AP.thetaLo=                    [4 6];
AP.thetaHi=                    [6 12];
AP.beta=                       [15 30];
AP.gamma=                      [40 90];
AP.ripple=                     [120 180];

% ------ cutoff frequencies for butterworth bandpass filters
AP.thetaCFreq=                 [4 12];
AP.thetaLoCFreq=               [4 6];
AP.thetaHiCFreq=               [6 12];
AP.betaCFreq=                  [15 30];
AP.gammaCFreq=                 [40 90];
AP.rippleCFreq=                [130 170];

% ------ steepness of filters in dB/octave 
% (see help in bafi.m for more detailed information)
AP.rs=                         50;

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% SPECTRAL & CORRELATION ANALYSIS PARAMETERS
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% number of points in elementary data segment on which both spectra and
% crosscorrelations are based.
% IMPORTANT: this is also the minimum acceptable duration of behavioral segments; i.e. if
% an animal performs any of the behaviors above for less than this time in a row, or
% two artifacts (rather, their pre- and post artifact windows) embrace an interval less
% than this time that behavioral segment is lost for analysis. 
% NOTES: 
% - must be a power of 2 
% - spectral resolution is sampling frequency/points per segment, e.g. 1000 Hz/2048 pts
%   would give a resolution of ~0.5 Hz
AP.ppSeg=                     4096;

% ------ Discrete Fourier Transform - DFT
% overlap between segments - anywhere between a fourth and half of the above will do the trick
AP.dftOlapPts=                 round(AP.ppSeg/3);
% this is the window by which each data segment will be multiplied prior to spectral
% analysis. 
AP.dftWin=                     'Hamming';

% ------ Cross correlations - CC
% lag of CC expressed in number of points 
AP.ccLagPts=                    1000;   
% scaling of CC, see help for function xxcorr
AP.ccScaleOpt=                 'coeff_ub';
% the interval in which to detect & save peaks in all CC obtained with theta signals
AP.thccw=                      [-500 500];
% the interval in which to detect & save peaks in the gamma CC 
AP.gaccw=                      [-100 100];
% CC theta,gammaEnv: the interval within which to pick the peaks in the AVERAGE cc.
% These peaks are the 'attractors' for cc peaks from individual segments. This interval 
% must be completely contained in the above AP.thccw
AP.cca=                        [-150 150];
% the significance of the CC between theta and gamma envelope is evaluated by
% computing the CC between theta and phase-shuffled versions of the gamma envelope.
% The parameter below specifies on how many shuffled versions of gammaEnv this
% evaluation should be based. Set to zero to skip shuffling.
AP.ccNShuffle=                  50;   

% ----- spectral peak frequency detection
% .. is accomplished as follows:
% the power spectra (averaged over all individual segments) are
% 1. multiplied by a function which greatly attenuates values outside the frequency
% range of interest (more precisely, this function is the DFT of a FIR bandpass filter, 
% i.e. multiplication with it is equivalent to (and, for that matter, computationally 
% more efficient than) having filtered the data with that filter)
% 2. smoothed with a triangular window.
% In the resulting spectrum, the most prominent peak is determined. The parameters
% below specify the details of these different steps

% frequency range (Hz) within which to search for (theta) peak in power spectrum. 
% The FIR filter mentioned above will have corner frequencies approximately matching 
% these values
AP.peakFRng=                     [2.5 14];

% FIR bandpass filter order; must be < AP.ppSeg. The higher this value the sharper 
% the filter's rolloff
AP.FIRfo=                        400;

% the power spectra will be smoothed by convoluting them with a triangular window.
% the parameter below specifies the half width of that window in Hz
AP.specKrnlHWid=                 .24;

% if multiple peaks are found in the frequency range of interest, their 'prominence'
% is evaluated by computing their integral in the range [peak location +/-AP.peakHW]
% (in Hz). Also, the frequency band [peak location of *principal* channel +/-AP.peakHW] 
% will be used to estimate average coherence in the theta range 
AP.peakHW=                        .5;

% end of AP definition
