function rmouse_detstreampeak(rawCh,strmType)
% this is a behavior-independent type of analysis, so treat the results as
% streams, putting them in a separate file (instead of appending to results
% variable r)
% Also, detect peaks over the whole length of the streams (as present in
% the corresponding .i16 files)

global AP WP 

switch strmType
  case 'theta'
    fnList=strvcat(rawCh.thetaFn);
    % this should ensure smoothing of edges with minimal distortion of the
    % signal proper
    fiFreq=AP.thetaCFreq(2)+5;
    fiRs=30;
    sFn=WP.thetaPFn;    
  case 'gammaEnv'
    fnList=strvcat(rawCh.gammaEnvFn);
    % the gamma envelope is pretty ragged, so smooth it vehemently
    fiFreq=AP.thetaCFreq(2);
    fiRs=40;    
    sFn=WP.gammaEnvPFn;
  otherwise
    error('illegal stream type');
end
    

% preallocate results vars: all LFP channels
posPA=cell(AP.nAllLFPCh,1);
posPT=posPA;
negPA=posPA;
negPT=posPA;
nPosP=zeros(AP.nAllLFPCh,1);
nNegP=zeros(AP.nAllLFPCh,1);

% loop over LFP channels to be analyzed
for chInd=1:AP.nLFPCh
  % d=strmread([AP.strmDir '\' rawCh(AP.LFPInd(chInd)).   theta    Fn],'nPts',inf,'verbose',0);
  d=strmread([AP.strmDir '\' fnList(chInd,:)],'nPts',inf,'verbose',0);  
  % it's a drag, but inevitable: re-filter the streams because channels
  % with a very low amplitude suffer from the low resolution of the int16,
  % causing evdeal to detect peaks where there are none in the smooth
  % version
  d=lofi(d,WP.osi,fiFreq,'rs',fiRs);
  n1=length(d);
  tmpr=evdeal(d,WP.osi,'allpeaks');
  posPA(AP.LFPccInd(chInd))=tmpr.posPeak;
  posPT(AP.LFPccInd(chInd))=tmpr.posPeakT;
  nPosP(AP.LFPccInd(chInd))=length(tmpr.posPeakT{:});
  negPA(AP.LFPccInd(chInd))=tmpr.negPeak;
  negPT(AP.LFPccInd(chInd))=tmpr.negPeakT;
  nNegP(AP.LFPccInd(chInd))=length(tmpr.negPeakT{:});
end
save(sFn,'posPA','posPT','nPosP','negPA','negPT','nNegP');

