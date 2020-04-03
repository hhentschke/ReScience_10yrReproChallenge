% This code looks at the phase lags introduced by the Neuralynx hardware
% filters

% sampling freq 1000 Hz, filter (Neuralynx) 1-325 Hz, gain 1000
% IN 0 is direct input into axon digitizer, IN 1 via Neuralynx system

% 0000  10 Hz
% 0001  50 Hz
% 0002  100 Hz
% 0003  ~2 Hz
% 0004  ~5 Hz
% 0005  8 Hz

[d,si]=abfload('2005-06-07 0005.abf');

d(:,2)=d(:,2)*-1;
cc(d(:,1),d(:,2),50);
grid on