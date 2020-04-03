% quick-and dirty data plotter and/or wav-player

% pfac>1 => play faster
pfac=10;

f='e:\rmouse\wt2141_02510\02510000.abf',
chans={'IN 5','IN 6','IN 7','IN 8','IN 9','IN 10','IN 11','IN 12','IN 13','IN 14','IN 15',}
chans={'IN 0','IN 1','IN 2','IN 3','IN 4','IN 5','IN 6','IN 7','IN 8','IN 9','IN 10','IN 11','IN 12','IN 13','IN 14','IN 15',}

% wild field spikes @ ca 20 Hz appear on chans IN 7 - IN 15
intv=[77 83];
% detail
intv=[79.5 79.7];

% another episode, with clearer lags between the spikes
intv=[193.7 195.7];
% detail
intv=[193.7 194.0]+.9;

% one episode with field spikes in dorsal electrodes
intv=[174.6 176.0];

[d,si]=abfload(f,'start',intv(1),'stop',intv(2),'channels',chans);
% ** invert data?

clf, orient tall

%d=hifi(d,si,5);

pllplot(d,'si',si,'spacing','fixed','dy',2.0);

title([f ', chans IN 0 - IN 15, intv=' num2str(intv,4) ' s'])
%print -djpeg90 p

%wavplay(d,1e6/si*pfac,'sync');
