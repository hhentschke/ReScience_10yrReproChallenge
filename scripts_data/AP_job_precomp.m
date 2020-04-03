% this is a partial AP (Analysis Parameter) set to be used in batch processing 

AP.job={...
  'det_artifacts',...       % detect artifacts in neural signals & generate time stamp list
  'gen_btsl',...            % generate behavioral time stamp list
  'chechao',...             % check proper order of channels by crosscorrelating raw data
  'filter&hilbert',...      % generate the various streams (theta, gamma, etc.)
  'thetaPeaks',...          % compute amplitude & time of occurrence of peaks in theta stream
  'seg_thetaPeakReg',...    % compute theta regularity based on variance of peak amplitude and inter-peak interval
  'specgram',...            % compute & plot spectrogram (time course of spectral power)
  'seg_rawSpecPow',...      % compute spectra & derived parameters
  'seg_gammaEnvSpecPow',... % compute spectra & derived parameters for gamma envelope (for each channel)  
  'seg_deltaCC',...         % compute delta crosscorrelations (all combinations of channels)
  'seg_thetaCC',...         % compute theta crosscorrelations (all combinations of channels)
  'seg_thetaLoEnvCC',...    % compute theta_lo envelope crosscorrelations (all combinations of channels)
  'seg_thetaHiEnvCC',...    % compute theta_hi envelope crosscorrelations (all combinations of channels)
  'seg_gammaCC',...         % compute gamma crosscorrelations (all combinations of channels)
  'seg_gammaEnvCC',...      % compute gamma envelope crosscorrelations (all combinations of channels)
  'seg_thetaGammaEnvCC',... % compute crosscorrelations between theta and gamma envelope (for each channel)
  'seg_deltaThetaEnvCC',... % compute crosscorrelations between delta and theta envelope (for each channel)
  'seg_cc_ms2rad',...       % convert theta-related phase lags from ms to radians, creating new results fields (NOTE: will be called automatically if any of 'seg_rawSpecPow','seg_thetaCC','seg_gammaEnvCC','seg_thetaGammaEnvCC' is called)
  'delStrms',...            % delete the large files for theta, gamma, etc streams
  'sumFig',...              % produce summary plots
  'rawPlot',...             % plot of raw data excerpts
  'csdPlot',...             % contour plot of current source densities from theta stream excerpts  
  'tcFig',...               % plots time course of segment-wise computed, selected parameters 
  'tcCorr',...              % correlates segment-wise computed, selected parameters with each other 
  'tcPrincComp',...         % computes & plots principal components from tupels of segment-wise computed, selected parameters 
  'tcThetaPeak',...         % spiel: produce blob-plots depicting time course of theta peaks
  'gen_audio',...           % spiel (& under construction): generate wav files encoding theta, gamma and ripple signals as drumbeats, pink noise and chirp, respectively
  'spiel'...                % playground
};
  
AP.job={...
  'det_artifacts',...       % detect artifacts in neural signals & generate time stamp list
  'gen_btsl',...            % generate behavioral time stamp list
  'chechao',...             % check proper order of channels by crosscorrelating raw data
  'specgram',...            % compute & plot spectrogram (time course of spectral power)
  'seg_rawSpecPow',...      % compute spectra & derived parameters
  'sumFig',...              % produce summary plots
};


% determines how to deal with previously computed results of the 'seg_' jobs:
% 'replace' means they will be discarded entirely and replaced by the results 
%   of the current run. This should be the default.
% 'append' means that the results of the current run will be appended to or 
%   overwrite results of previous runs. 
% *** see NOTE of AP.job above! ***
AP.saveMode=                'replace';'append';


% save plot(s) to file in what format? '-dpsc2' (color postscript level 2) is
% strongly recommended. multiple formats can be specified; in this case they must 
% be listed in a cell array, like {'-dpsc2','-djpeg90'}. set to [] if plots 
% shall not be saved.
AP.printas=                  {'-djpeg90'};[];


% the significance of the CC between theta and gamma envelope is evaluated by
% computing the CC between theta and phase-shuffled versions of the gamma envelope.
% The parameter below specifies on how many shuffled versions of gammaEnv this
% evaluation should be based. Set to zero to skip shuffling.
AP.ccNShuffle=                  0;   
