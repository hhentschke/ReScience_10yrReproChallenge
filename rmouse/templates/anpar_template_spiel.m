% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%  PLAYGROUND: AUDIO 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% channels from which to generate audio output
AP.rawChAudioNm=             {'Lynx5','Lynx9'};

% length of segment
AP.rawAudioExcerpt=          'full'; [0 10];

% sampling frequency (Hz) of audio stream. Realize that 
% - neuronal data will be upsampled to this frequency, resulting in a stream
%   several times as large as the corresponding portion of neuronal data
% - the basic wav file(s) below will be played with that frequency, irrespective of 
%   whether it corresponds to its inherent one or not
AP.AudioFs=                   11025;

% name of the file with the basic wave
AP.wffn='d:\p&p\sounds\gen\drum1_mod.wav';

% that wave file's 'acoustic center' in points
AP.wCenter=                     82;

% the relative magnitudes of theta, gamma, ripple sounds (in this order) in wav file
AP.wMagScale=                 [1 .2 .6];

% modulation depth in 'arbitrary waveform frequency modulation' mode
% (a value between 0 and 1): amplitudes of field potentials will be 
% mapped to frequencies such that peaks correspond to high frequencies and troughs
% to low frequencies. The closer AP.modDepth is to 1, the more pronounced this
% modulation will be
AP.modDepth=                   .7;
