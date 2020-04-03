% ************************************************************************************
% This AP contains parameter values that should be common to every data set 
% processed with shorter segments for time course plots/analyses
% ************************************************************************************

global AP

AP.resFn=                    [DS.abfFn '_proc_iStateGamma'];
AP.saveMode=                'replace'; 'append';
AP.logFn=                    [DS.abfFn '_log_iStateGamma.txt'];
AP.printas=                  {'-djpeg90'};[];

% ms, the trigger signal's rise time from base line to threshold (in behavioral 
% scoring data) ==>> make 100% sure this is the slowest rise time the trigger pulses 
% will ever show, otherwise categorization of segments will get messed up seriously
AP.trigTau=                  20;

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% FREQUENCIES & FILTERS
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% ------ definition of frequency bands (Hz)
AP.delta=                      [1 4];
AP.theta=                      [5 12];
AP.thetaLo=                    [4 6];
AP.thetaHi=                    [6 12];
AP.beta=                       [15 30];
AP.gamma=                      [40 90];
AP.ripple=                     [120 180];

% ------ cutoff frequencies for butterworth bandpass filters
AP.deltaCFreq=                 [1 3.5];
AP.thetaCFreq=                 [5 12];
AP.thetaLoCFreq=               [4 6];
AP.thetaHiCFreq=               [6 12];
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
AP.ppSeg=                     64;

% ------ Discrete Fourier Transform - DFT
% overlap between segments - anywhere between a fourth and half of the above will do the trick
AP.dftOlapPts=                 round(AP.ppSeg/2);
AP.dftOlapPts=                 54;
% this is the window by which each data segment will be multiplied prior to spectral
% analysis. 
AP.dftWin=                     'hamming';

% ------ Cross correlations - CC
% lag of CC expressed in number of points 
AP.ccLagPts=                    45;   
% scaling of CC, see help for function xxcorr
AP.ccScaleOpt=                 'coeff';'coeff_ub';
% the interval in which to detect & save peaks in all CC obtained with delta signals
AP.deccw=                      [-395 395];
% the interval in which to detect & save peaks in all CC obtained with theta signals
AP.thccw=                      [-300 300];
% the interval in which to detect & save peaks in the gamma CC 
AP.gaccw=                      [-40 40];
% CC theta,gammaEnv: the interval within which to pick the peaks in the AVERAGE cc.
% These peaks are the 'attractors' for cc peaks from individual segments. This interval 
% must be completely contained in the above AP.thccw
AP.cca=                        [-150 150];

% the significance of the CC between theta and gamma envelope is evaluated by
% computing the CC between theta and phase-shuffled versions of the gamma envelope.
% The parameter below specifies on how many shuffled versions of gammaEnv this
% evaluation should be based. Set to zero to skip shuffling.
AP.ccNShuffle=                  0;   

