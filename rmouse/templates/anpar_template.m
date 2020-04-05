% ************************************************************************************
%                  This AP is the master template for other APs                 
% ************************************************************************************
% - it contains ALL fields of AP
% - it contains the most comprehensive and up-to-date information about the fields
% - do not mess around with it - rmouse_APcheck calls it to compare fields
% - the values of the individual fields are suggestions/have proved to work but are 
%   overridden by data set-specific values when rmouse is called
% ************************************************************************************


global AP DS

% invoke the data set-generating routine
% dset_template;

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%               JOB 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% The partial analysis job(s) to be performed. Further information on how they work:
% - code in rmouse.m and m-files called by rmouse.m
% - manual for rmouse
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
  'tcParCorr',...           % computes & plots correlations between segment-wise values of different results
  'tcThetaPeak',...         % produce blob-plots depicting time course of theta peaks (the plots have an aesthetic appeal and may reveal variability of theta)
};

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%  LENGTH OF DATA; CHANNELS; FILE NAMES; PRINT ETC
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% time block(s) to be processed in s (from t=0=start of recording).
% May be given as 'full'
AP.rawExcerpt=              'full';[0 200];

% the channels to be analyzed (must be a subset of DS.rawChan)
% NOTE 1: non-neuronal channels (scoring, wheel velocity, etc.) must be excluded
% NOTE 2: specify full names as given in the abf file header in cell 
% array of strings, like {'IN 0','IN 10'}
AP.rawChAnNm=               DS.rawCh(:,1);

% the 'principal' channel (serving as a reference for e.g. CC computations) 
% Same rules as above apply
AP.rawChPrincNm=             {'IN 7'};

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

% directory for streams files - will be generated automatically 
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
AP.bScoreFn=                '02313000a_bscore';

% name of the behavioral scoring channel in abf file (will be ignored if text file 
% is used)
AP.rawChScoreNm=            {'<enter scoring channel name>'};

% s, point in time where behavioral recording was started, relative to the beginning
% of the neural recording (as indicated by timer in camera view). positive values
% mean the behavioral recording started after the neural recording.
AP.bbor=                     0;

% ms, the assumed average delay between the onset of a behavior and the experimenter's 
% reaction (pressing a button)
AP.bsRt=                       300;

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
    'grooming',  [1.0 2.4],  [.5 .4 1],  '^';...
    'exploring', [6   9],  [.5 1 .5],   'o';...
    'immobile',  [2.5 4],  [.6 .4 .1],  's';...
    'bad',       [-1 -9],  [.9 0 0],    '+'...
};

% ms, the trigger signal's rise time from base line to threshold (in behavioral 
% scoring data) ==>> make 100% sure this is the slowest rise time the trigger pulses 
% will ever show, otherwise categorization of segments will get messed up seriously
AP.trigTau=                  20;

% Hz, corner frequency of highpass filter to apply to behavioral scoring data 
% (in some instances the baseline drifts slowly). Set to [] to skip filtering.
AP.bsCFreq=                  .1;

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%        DETECTION OF ARTIFACTS
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% the neuronal channel(s) to be monitored for periods of bad signal quality (noise) 
% and/or other unwanted periods. Must be a subset of AP.rawChAnNm
AP.rawChMonNm=               {'IN 8'};

% Hz, corner frequencies of lowpass filter to apply to channel monitored for artifacts
% set to [] if no filter shall be applied
AP.afCFreq=                   [];

% mV, threshold(s) for detection of artifacts on raw channels specified in 
% AP.rawChMonNm 
AP.afThresh=                  [-1.0];

% seconds, the interval around an artifactual event/period (noise) to omit from analysis
AP.afWin=                      [-.5 .5];

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% FREQUENCIES & FILTERS
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ------ list of streams to be generated (=written into binary files)
% ** Note: all of the envelope streams require the 'ancestor' stream to be
% generated as well, and must be listed right below it)
AP.strm={...
  'sansDelta',...
  'delta',...
  'theta',...
  'thetaEnv',...
  'thetaLo',...
  'thetaLoEnv',...
  'thetaHi',...
  'thetaHiEnv',...
  'beta',...
  'gamma',...
  'gammaEnv',...
  'gammaNarrow',...
  'gammaNarrowEnv',...
  'ripple',...
  'rippleEnv'};

% ------ definition of frequency bands (Hz)
AP.sansDelta=                  [4 400];
AP.delta=                      [1 4];
AP.theta=                      [4 12];
AP.thetaLo=                    [4 6];
AP.thetaHi=                    [6 12];
AP.beta=                       [15 30];
AP.gamma=                      [40 90];
AP.gammaNarrow=                [50 80];
AP.ripple=                     [120 180];

% ------ cutoff frequencies for butterworth bandpass filters
AP.sansDeltaCFreq=             [4 400];
AP.deltaCFreq=                 [1 3.5];
AP.thetaCFreq=                 [4 12];
AP.thetaLoCFreq=               [4 6];
AP.thetaHiCFreq=               [6 12];
AP.betaCFreq=                  [15 30];
AP.gammaCFreq=                 [40 90];
AP.gammaNarrowCFreq=           [50 80];
AP.rippleCFreq=                [130 170];


% ------ center frequency for butterworth bandstop filter. If nonempty, the gamma stream 
% will be filtered at corner frequencies +/-1 this value with high rolloff (set to [] 
% to prevent bandstop filtering)
AP.lineFreq=                   [60];

% ------ steepness of filters in dB/octave 
% (see help in bafi.m for more detailed information)
AP.rs=                         40;

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% SPECTRAL & CORRELATION ANALYSIS PARAMETERS
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% number of points in elementary data segment on which frequency spectra,
% and cross correlations and many other results are based.
% IMPORTANT: this is also the minimum acceptable duration of behavioral
% segments; i.e. if an animal shows any behavior for less than this time
% interval, or two artifacts (rather, their pre- and post artifact windows)
% embrace an interval less than AP.ppSeg, that behavioral segment is lost
% for analysis. 
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
AP.dftWin=                     'hamming';

% ------ Cross correlations - CC
% lag of CC expressed in NUMBER OF POINTS 
AP.ccLagPts=                    1000;   
% scaling of CC, see help for function xxcorr
AP.ccScaleOpt=                 'coeff';
%           **** the following parameters are specified in ms ****
% the interval in which to detect & save peaks in all CC obtained with delta signals
AP.deccw=                      [-500 500];
% the interval in which to detect & save peaks in all CC obtained with theta signals
AP.thccw=                      [-300 300];
% the interval in which to detect & save peaks in the gamma CC 
AP.gaccw=                      [-35 35];
% the interval in which to detect & save peaks in the NARROW gamma CC 
AP.gaNaccw=                    [-25 25];
% CC theta,gammaEnv: the interval within which to pick the peaks in the AVERAGE cc.
% These peaks are the 'attractors' for cc peaks from individual segments. This interval 
% must be completely contained in the above AP.thccw
AP.cca=                        [-150 150];
% the significance of the CC between theta and gamma envelope is evaluated by
% computing the CC between theta and phase-shuffled versions of the gamma envelope.
% The parameter below specifies on how many shuffled versions of gammaEnv this
% evaluation should be based. Set to zero to skip shuffling.
AP.ccNShuffle=                  0;   
% there are a few cases where the detection of prominent peaks in the
% theta-gammaEnv CC gets 'derailed' (lags of detected peaks are off),
% mostly due to either or a combination of very small peaks or
% inter-electrode phase jumps > half a theta period (in many cases this is
% caused by excessive line hum. Try getting rid of it via bandstop
% filtering (see AP.lineFreq above)). In these cases, you can specify, for
% each channel, the 'attractor' lag. Thus, in each crosscorrelation
% segment, the peak closest to this lag will be picked (in 'normal' mode,
% the attractor lag corresponds to the lag of the peak which is closest to
% the attractor lag of a neighboring electrode. There is always only one
% neighboring electrode for which the attractor lag had been calculated,
% and the whole procedure starts either with the first channel or a channel
% of your choice (see AP.OPT_thgaeCCStartChNm below). This parameter is
% OPTIONAL (i.e. it need not be specified). If specified, the number of
% lags must correspond to the number of ALL LFP channels TO BE ANALYZED
% (AP.rawChAnNm).
AP.OPT_thgaeCCAttractLag=[50 50 50 30 -20 -20 0 0 10 60 60 60 60 60 60];
% the channel to initiate computation of theta-gammaEnv CC with (by
% default, it is the first channel). This parameter is OPTIONAL (i.e. it
% need not be specified). 
AP.OPT_thgaeCCStartChNm=        {'IN 7'};

% same parameters for deltaThetaEnv
AP.OPT_detheCCAttractLag=[50 50 50 30 -20 -20 0 0 10 60 60 60 60 60 60];
AP.OPT_detheCCStartChNm=        {'IN 7'};


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
AP.peakFRng=                     [5 12];

% approximate corner frequencies of FIR filter mentioned above 
AP.FIRcf=                        [2 12];

% FIR bandpass filter order; must be < AP.ppSeg. The higher this value the sharper 
% the filter's rolloff
AP.FIRfo=                        400;

% the power spectra will be smoothed by convoluting them with a triangular window.
% the parameter below specifies the half width of that window in Hz. With a segment 
% length of 4096 points and a sampling frequency of ~1000 Hz (which is the
% case in most recordings), a half-width of 0.24 amounts to a 3-point
% triangular window
AP.specKrnlHWid=                 .24;

% if multiple peaks are found in the frequency range of interest, their 'prominence'
% is evaluated by computing their integral in the range [peak location +/-AP.peakHW]
% (in Hz). 
AP.peakHW=                        .5;

% end of AP definition
% ----------------------------------------------------------------------------------------
