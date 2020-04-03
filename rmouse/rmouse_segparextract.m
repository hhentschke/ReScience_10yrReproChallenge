function [T,Y,B,ub,ubix]=rmouse_segparextract(rv,behav,dDist)
% ------- extracts selected set of segment-wise results parameters for
% i) plots of these parameters, ii) principal components analysis with the
% aim of determining the segregation between behavioral states

global AP WP r logstr

% improvements: 
% - check whether channel is in the set of those computed
% - avoid repetitive execution (check whether vars exist)


nrv=length(rv);
[doff,ccRefChInd]=min(abs(WP.elx(AP.LFPccInd)+dDist));
ccRefChIdx=AP.LFPIdx(ccRefChInd);
% §§§ index wrong: WP.elx contains values for all LFPs while both
% ccRefChIdx and AP.pcIdx refer to ALL channels - change this
% dDist=diff(WP.elx([ccRefChIdx AP.pcIdx]));

% load r only if it that is not done yet 
if isempty(r), 
  if exist([AP.resFn '.mat'],'file')
    load(AP.resFn); 
  else
    logstr{end+1}='results file does not exist';
    warning(logstr{end});
  end
end


% variables collecting all data points 
Y=[]; T=[]; B=[];
ub=[]; ubix={[]};
% loop over behaviors
for i=1:length(r)
  if ~isempty(r(i).iPts) & strmatch(r(i).segmentType(:,1),behav)
    % time axis - minutes, please
    t=discrete2cont(round(mean(r(i).iPts,2)),WP.osi*.001)/6e4;
    nt=length(t);
    T=[T;t];
    % behavior code = index to r
    B=[B; repmat(i,nt,1)];
    ub=[ub i];
    if isempty(ubix{1})
      ubix{1}=1:nt;
    else
      ubix{end+1}=[1:nt]+ubix{end}(end);
    end
    Y=[Y; repmat(nan,nt,nrv)];
    for rvi=1:nrv
      eval(['vex=isfield(r(i),''' rv{rvi} ''');']);
      if vex, eval(['y=r(i).' rv{rvi} ';']);
      else y=[];
      end
      % select extraction method depending on results var chosen 
      switch rv{rvi}
        % the within-site vars (strictly speaking the CVs are not
        % within-site vars because they are not the result of computations
        % among different streams from one site, but they are stored in the
        % same way as the former)
        case {'thgaeCCPeak','thgaeCCPeakT','thgaeCCPeakPhase',...
            'thgaeCCPosPeakDecay','thgaeCCZScore','thgaeCCZTestP',...
            'thPosPeakCvA','thPosPeakCvIPI','thNegPeakCvA','thNegPeakCvIPI',...
            'gaePosPeakCvA','gaePosPeakCvIPI'}
          if ~isempty(y)
            y=y(:,AP.LFPpcInd1);
          else
            y=[];
          end
        % autos 
        case {'rawDePE','rawThPE','rawThNarrowPE','rawBePE','rawGaPE',...
            'rawGaNarrowPEMn','rawRiPE',...
            'rawGaCentroid',...
            'rawPPeak','rawPPeakT',...
            'thCCPosPeakDecay',...
            'gaePPeak','gaePPeakT',...
            'gaeThPEMn','gaeThNarrowPEMn'}
          if length(diag(y))>=1
            y=y{AP.LFPpcInd2,AP.LFPpcInd2};
          else
            y=[];
          end
        % all whose cross-measures shall be plotted
        otherwise
          if length(diag(y))>=1
            % y=y{ccRefChIdx,AP.LFPpcInd2};
            y=y{ccRefChInd,AP.LFPpcInd2};
          else
            y=[];
          end
      end % switch
      if isempty(y)
        disp(['field ' rv{rvi} ' does not exist or contains no data']);
      else
        % collect values
        Y(end-nt+1:end,rvi)=makecol(y);
      end
    end % for:rv
  end % if:~isempty r.ipts & behav requested
end % for:behaviors
