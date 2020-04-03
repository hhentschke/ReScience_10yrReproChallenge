function rmouse_cc_ms2rad
% This piece of code converts phase lags of theta CC, gammaEnv CC and 
% theta-gammaEnv CC from ms to radians, generating new fields of r in the process.
% The inverse of peak theta frequency of the principal channel (as obtained 
% from spectral analysis) is regarded as one period

global DS AP WP r 


% due to a bug in rmouse_ccintrasite means & stdev of parameters below had
% to be recomputed
% for ggg=3:4
%   r(ggg).thgaeCCPeakMn=nanmean(r(ggg).thgaeCCPeak);
%   r(ggg).thgaeCCPeakStd=nanstd(r(ggg).thgaeCCPeak);
%   r(ggg).thgaeCCPeakTMn=nanmean(r(ggg).thgaeCCPeakT);
%   r(ggg).thgaeCCPeakTStd=nanstd(r(ggg).thgaeCCPeakT);  
% end


% check whether peak freq had been determined before
if isfield(r,'rawPMnPeakT')
  peakFindAttempt=1;
else
  disp(['** warning (' mfilename '): theta-related phase lags cannot be converted to rad because results from spectral analysis are missing']);
  peakFindAttempt=0;
  cFac=nan;
end
% check whether fields to be converted exist
fexist(1)=isfield(r,'thCCPeakTMn');
fexist(2)=isfield(r,'gaeCCPeakTMn');
fexist(3)=isfield(r,'thgaeCCPeakTMn');
fexist(4)=isfield(r,'thgaeCCEnvPeakTMn');

for i=1:length(r)
  if ~isempty(r(i).iPts)
    if peakFindAttempt
      cFac=2*pi*0.001*r(i).rawPMnPeakT{AP.pcInd,AP.pcInd};
      cFacOK=1;
      if ~isfinite(cFac)
        % If no peak was found (cFac=nan) we have a pathological case (psd
        % with a slope < 0 over whole theta range)
        cFacOK=0;
        disp(['** warning (' AP.segmentType{i,1} '): theta-related phase lags cannot be converted to rad because no theta peak was found']);
      end
    else
      cFacOK=0;
    end
    % --- theta CC
    if fexist(1)
      r(i).thCCPeakPhaseMn=WP.ccDerTemplate;
      r(i).thCCPeakPhaseStd=WP.ccDerTemplate;
      if cFacOK
        for g=AP.trixie'
          r(i).thCCPeakPhaseMn{g}=r(i).thCCPeakTMn{g}*cFac;
          r(i).thCCPeakPhaseStd{g}=r(i).thCCPeakTStd{g}*cFac;
        end
      end
    else
      disp([AP.segmentType{i,1} ': theta phase lags not converted to rad because crosscorrelation results do not exist']);
    end
    % --- gammaEnv CC
    if fexist(2)
      r(i).gaeCCPeakPhaseMn=WP.ccDerTemplate;
      r(i).gaeCCPeakPhaseStd=WP.ccDerTemplate;
      if cFacOK
        for g=AP.trixie'
          r(i).gaeCCPeakPhaseMn{g}=r(i).gaeCCPeakTMn{g}*cFac;
          r(i).gaeCCPeakPhaseStd{g}=r(i).gaeCCPeakTStd{g}*cFac;
        end
      end
    else
      disp([AP.segmentType{i,1} ': gammaEnv phase lags not converted to rad because crosscorrelation results do not exist']);
    end
    if fexist(3)
      % --- with theta-gammaEnv things are much simpler because the results are in an array
      r(i).thgaeCCPeakPhaseMn=r(i).thgaeCCPeakTMn*cFac;
      r(i).thgaeCCPeakPhaseStd=r(i).thgaeCCPeakTStd*cFac;
    else
      disp([AP.segmentType{i,1} ': theta-gammaEnv phase lags not converted to rad because crosscorrelation results do not exist']);
    end
    if fexist(4)
      % --- envelope of theta-gammaEnv cc
      r(i).thgaeCCEnvPeakPhaseMn=r(i).thgaeCCEnvPeakTMn*cFac;
      r(i).thgaeCCEnvPeakPhaseStd=r(i).thgaeCCEnvPeakTStd*cFac;
    else
      disp([AP.segmentType{i,1} ': theta-gammaEnv (env) phase lags not converted to rad because crosscorrelation results do not exist']);
    end
  end
end