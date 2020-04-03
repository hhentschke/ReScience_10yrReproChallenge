global AP

% invoke the data set-generating routine
dset_2004_07_30_0006;

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%               JOB 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% the partial analysis job(s) to be performed. See the analysis manual for an
% explanation of the workings behind it. 
% A few very basic notes/reminders:
% - when running a subset of the partial jobs, make sure the basic ones,
%   upon which the others depend, are included:
% -- all jobs starting with 'seg_' do computations on (behaviorally scored) 
%    segments of data and therefore must always be preceded by job 'gen_btsl'
% -- 'filter&hilbert' generates the numerous intermediate files containing
%     gamma, theta etc. streams and needs to be called only once (unless the
%     filter settings change or the streams were deleted)
% - do not change order of jobs
% - when you change parameters from one run to the next, run the whole
%   analysis from scratch, even if that amounts to a lot of redundant computations,
%   unless you know exactly which jobs are unaffected by the changed parameter(s)

AP.job={...
  'det_artifacts',...       % detect artifacts in neural signals & generate time stamp list
  'gen_btsl',...            % generate behavioral time stamp list
  'chechao',...             % check proper order of channels by crosscorrelating raw data
  'filter&hilbert',...      % generate the various streams (theta, gamma, etc.)
  'thetaPeaks',...          % compute amplitude & time of occurrence of peaks in theta stream
  'seg_rawSpecPow',...      % compute spectra & derived parameters
  'seg_thetaCC',...         % compute theta crosscorrelations (all combinations of channels)
  'seg_gammaCC',...         % compute gamma crosscorrelations (all combinations of channels)
  'seg_thetaEnvCC',...      % compute theta envelope crosscorrelations (all combinations of channels)
  'seg_gammaEnvCC',...      % compute gamma envelope crosscorrelations (all combinations of channels)
  'seg_thetaGammaEnvCC',... % compute crosscorrelations between theta and gamma envelope (for each channel)
  'delStrms',...            % delete the large files for theta, gamma, etc streams
  'sumFig',...              % produce summary plots
  'rawPlot',...             % plot of raw data excerpts
  'csdPlot',...             % contour plot of current source densities from theta stream excerpts  
  'tcFig',...               % plots time course of segment-wise computed, selected parameters 
  'tcPrincComp',...         % computes & plots principal components from tupels of segment-wise computed, selected parameters 
};

AP.job={...
  'det_artifacts',...       % detect artifacts in neural signals & generate time stamp list
  'gen_btsl',...            % generate behavioral time stamp list
  'seg_gammaCC',...         % compute gamma crosscorrelations (all combinations of channels)
  'sumFig',...              % produce summary plots
};

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%  LENGTH OF DATA; CHANNELS; FILE NAMES; PRINT ETC
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% time block(s) to be processed in s (from t=0=start of recording).
% May be given as 'full'
AP.rawExcerpt=              'full';[0 200];

% the channels to be analyzed (must be a subset of DS.rawChan)
% NOTE: specify full names as given in the abf file header in cell 
% array of strings, like {'IN 0','IN 10'}
AP.rawChAnNm=               DS.rawCh(1:15,1);

% the 'principal' channel (serving as a reference for e.g. CC computations) 
% Same rules as above apply
AP.rawChPrincNm=             {'IN 8'};

% determines how to deal with previously computed results of the 'seg_' jobs:
% 'replace' means they will be discarded entirely and replaced by the results 
%   of the current run. This should be the default.
% 'append' means that the results of the current run will be appended to or 
%   overwrite results of previous runs. 
% *** see NOTE of AP.job above! ***
AP.saveMode=                'append';'replace';

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

% save plot(s) to file in what format? '-dpsc2' (color postscript level 2) is
% strongly recommended for high-quality figures, whereas '-djpeg90' (jpeg with 
% quality level of 90) is more efficient for everyday use. Multiple formats can 
% be specified; in this case they must be listed in a cell array, like 
% {'-dpsc2','-djpeg90'}. Set to [] if plots shall not be saved.
AP.printas=                  {'-djpeg90'};

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
AP.bScoreFn=                '2004-07-30 0006';

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
AP.rawChMonNm=               {'IN 3'};

% mV, threshold(s) for detection of artifacts on raw channels specified in 
% AP.rawChMonNm 
AP.afThresh=                  [2];

% Hz, corner frequencies of lowpass filter to apply to channel monitored for artifacts
% set to [] if no filter shall be applied
AP.afCFreq=                   [];

% seconds, the interval around an artifactual event/period (noise) to omit from analysis
AP.afWin=                      [-2 2];

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% FREQUENCIES & FILTERS
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ------ definition of frequency bands (Hz)
AP.delta=                      [1 4];
AP.theta=                      [4 12];
AP.gamma=                      [30 90];
AP.ripple=                     [120 180];

% ------ cutoff frequencies for butterworth bandpass filters
AP.deltaCFreq=                 [1 4];
AP.thetaCFreq=                 [4 12];
AP.gammaCFreq=                 [30 90];
AP.rippleCFreq=                [130 170];

% ------ center frequency for butterworth bandstop filter. If nonempty, the gamma stream 
% will be filtered at corner frequencies +/-1 this value with high rolloff (set to [] 
% to prevent bandstop filtering)
AP.lineFreq=                   [60];

% ------ steepness of filters in dB/octave 
% (see help in bafi.m for more detailed information)
AP.rs=                         50;

% % ------ parameter files for filters - not yet implemented
% % NOTE: this option overrides the edge frequencies specified above!
% % Set to {} to deactivate
% AP.thetaFparFn=                {}; %{'d:\filt\theta_a.asc','d:\filt\theta_b.asc'};
% AP.gammaFparFn=                {}; %{'d:\filt\gamma_a.asc','d:\filt\gamma_b.asc'};

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
AP.ccNShuffle=                  0;   

% there are a few cases where the detection of prominent peaks in the theta-gammaEnv
% CC gets 'derailed' (lags of detected peaks are off), mostly when very small peaks 
% or inter-electrode phase jumps > half a theta period occur. In these cases, you can 
% specify, for each channel, the 'attractor' lag. thus, in each crosscorrelation segment, 
% the peak closest to this lag will be picked (in 'normal' mode, the attractor lag 
% corresponds to the lag of the peak which is closest to the attractor lag of a neighboring 
% electrode. There is always only one neighboring electrode for which the attractor lag had 
% been calculated, and the whole procedure starts with a channel with the strongest average 
% cross spectral power in theta and gamma frequency bands).
% This parameter is OPTIONAL (i.e. it need not be specified). If specified, the number of
% lags must correspond to the number of ALL LFP channels.
AP.OPT_thgaeCCAttractLag=[50 50 50 50 40 -20 -20 zeros(1,8)];

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
AP.peakFRng=                     [4 15];

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
