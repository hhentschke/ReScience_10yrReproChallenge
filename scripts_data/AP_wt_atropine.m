% ************************************************************************************
% This AP contains parameter values that should be common to every data set
% ************************************************************************************

global AP DS


% ms, the trigger signal's rise time from base line to threshold (in behavioral 
% scoring data) ==>> make 100% sure this is the slowest rise time the trigger pulses 
% will ever show, otherwise categorization of segments will get messed up seriously
AP.trigTau=                  20;

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% FREQUENCIES & FILTERS
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ------ list of streams to be generated (=written into binary files)
% ** Note: all of the envelope streams require the 'ancestor' stream to be
% generated as well, and must be listed right below it)
AP.strm={...
  'sansDelta',...
  'theta',...
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
AP.thetaHi=                    [8 11];
AP.beta=                       [15 30];
AP.gamma=                      [30 90];
AP.gammaNarrow=                [65 90];
AP.ripple=                     [120 180];

% ------ cutoff frequencies for butterworth bandpass filters
AP.sansDeltaCFreq=             [4 400];
AP.deltaCFreq=                 [1 3.5];
AP.thetaCFreq=                 [4 12];
AP.thetaLoCFreq=               [4 6];
AP.thetaHiCFreq=               [8 11];
AP.betaCFreq=                  [15 30];
AP.gammaCFreq=                 [30 90];
AP.gammaNarrowCFreq=           [65 90];
AP.rippleCFreq=                [130 170];

% ------ steepness of filters in dB/octave 
% (see help in bafi.m for more detailed information)
AP.rs=                         40;


% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% SPECTRAL & CORRELATION ANALYSIS PARAMETERS
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
AP.ppSeg=                     4096;

% ------ Discrete Fourier Transform - DFT
% overlap between segments - anywhere between a fourth and half of the above will do the trick
AP.dftOlapPts=                 round(AP.ppSeg/3);
% this is the window by which each data segment will be multiplied prior to spectral
% analysis. 
AP.dftWin=                     'hamming';

% ------ Cross correlations - CC
% lag of CC expressed in number of points 
AP.ccLagPts=                    1000;   
% scaling of CC, see help for function xxcorr
AP.ccScaleOpt=                 'coeff_ub';
% the interval in which to detect & save peaks in all CC obtained with delta signals
AP.deccw=                      [-500 500];
% the interval in which to detect & save peaks in all CC obtained with theta signals
AP.thccw=                      [-300 300];
% the interval in which to detect & save peaks in the gamma CC 
AP.gaccw=                      [-35 35];
% the interval in which to detect & save peaks in the NARROW gamma CC 
AP.gaNaccw=                    [-25 25];
% triple purpose window for CC theta,gammaEnv: 
% a) for the principal channel, this is the interval within which to pick the most 
% prominent peak in the AVERAGE cc. This peak is the 'attractor' for cc peaks 
% from individual segments
% b) all other channels: principally the same, only the window is centered 
% around each channel's attractor peak
% c) all channels, INDIVIDUAL segments: window centered around attractor
% lag whithin which to search for peaks
AP.cca=                        [-100 100];

% ----- spectral peak frequency detection

% frequency range (Hz) within which to search for (theta) peak in power spectrum. 
AP.peakFRng=                     [5 12];

% approximate corner frequencies of FIR filter mentioned above 
AP.FIRcf=                        [2.5 13.5];

% FIR bandpass filter order; must be < AP.ppSeg. The higher this value the sharper 
% the filter's rolloff
AP.FIRfo=                        420;

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
