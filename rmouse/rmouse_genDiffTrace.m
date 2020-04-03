function rmouse_genDiffTrace
% ** function rmouse_genDiffTrace
% generates a time series to be used for the detection of epileptiform FP
% 'spikes' and artifacts:
% - filter
% - determine abs(hilbert(slope))
% - normalize
% - average
% The result will be written in abf-compatible mat-format (which can be
% read by matDload.m)

global DS AP WP

% all channels from alveus to HF 
locChIx=fliplr(AP.pcInd:-1:max(1,AP.pcInd-6));
% name of results file
strmFn=[AP.strmDir '\' DS.abfFn '_diff'];

if exist([DS.dpath '\' DS.abfFn '.mat'],'file')
  load([DS.dpath '\' DS.abfFn '.mat'],'abfi');
  [d,si]=matDload([DS.dpath '\' DS.abfFn '.mat'],'start',0,'stop','e','channels',AP.rawChAnNm(locChIx));

elseif exist([DS.dpath '\' DS.abfFn '.raw'],'file')
  rawFInfo=rawinfo('filename',[DS.dpath '\' DS.abfFn '.raw'],'print','no');
  d=rawload([DS.dpath '\' DS.abfFn '.raw'],AP.rawChAnNm(locChIx),'full',rawFInfo);
  % put into array and convert to mV
  d=cat(2,d{:})/1000;
  abfi=rawfi2genfi(rawFInfo,'channels',AP.rawChAnNm);
  si=abfi.si;
else
  [nix,nix2,abfi]=abfload([DS.dpath '\' DS.abfFn '.abf'],'info');
  % abfi=abfinfo([DS.dpath '\' DS.abfFn '.abf']);
  [d,si]=abfload([DS.dpath '\' DS.abfFn '.abf'],'start',0,'stop','e','channels',AP.rawChAnNm(locChIx));
end
[n1 n2]=size(d);
% filter: eliminate hi freqs > 40 Hz
d=lofi(d,si,40);
% envelope of slope: channel by channel to reduce memory consumption
for g=1:n2
  % unit: mV/ms
  d(:,g)=[0; abs(hilbert(diff(d(:,g))/(si/1000)))];
  % (append one data point at beginning so that number of points is equal to
  % raw data)
end
% normalize by 90th percentile
nrm=prctile(d,90,1);
d=d./repmat(nrm,n1,1);
% now average
d=single(mean(d,2));
% ** generate file information structure fi containing some essentials
% derived from abfi
fi.si=abfi.si;
fi.recChNames={'d'};
fi.dataPtsPerChan=abfi.dataPtsPerChan;
fi.dataPts=abfi.dataPtsPerChan;
fi.recTime=abfi.recTime;

save(strmFn,'d','fi');