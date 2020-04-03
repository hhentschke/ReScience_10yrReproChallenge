function rmouse_segoscillregularity(strmType,rawCh,STshort)
% ** function rmouse_segoscillregularity(strmType,rawCh,STshort)
% This routine computes coefficients of variation (cv) of peak amplitude
% and inter-peak-intervals (ipi) of an oscillatory input stream (theta or
% gamma envelope) on a segment by segment basis. The idea is to find out
% how regular the oscillation is.

% export functionality 'repaired' and extended Sep 2010

global AP WP r logstr

% set to true if behaviorally sorted peaks are to be exported (slows
% execution down)
doExport=false;
% minimally acceptable ipi
minIPI=5;

% load peaks (time is ms)
switch strmType
  case 'theta'
    pFn=WP.thetaPFn;
  case 'gammaEnv'
    pFn=WP.gammaEnvPFn;
  otherwise
    error('illegal stream type');
end
load(pFn,'negPA','negPT','posPA','posPT');

% --- first loop: pre-process peak data independent of behavior
negPIPI=cell(size(negPA));
posPIPI=negPIPI;
logstr{end+1}='preprocessing..';
disp(logstr{end});
for chaI=1:AP.nLFPCh
  logstr{end+1}=[rawCh(AP.LFPInd(chaI)).nm ' ..'];
  disp(logstr{end});
  % assign to local variables
  npA=negPA{AP.LFPccInd(chaI)};
  npT=negPT{AP.LFPccInd(chaI)};
  ppA=posPA{AP.LFPccInd(chaI)};
  ppT=posPT{AP.LFPccInd(chaI)};

  % kick out peaks outside bounds defined by AP.rawExcerpt (speeds up
  % computations for small excerpts)
  npIx=npT>=WP.boe & npT<WP.eoe;
  ppIx=ppT>=WP.boe & ppT<WP.eoe;
  logstr{end+1}=[int2str(length(find(npIx))) ' negative peaks out of ' int2str(length(npIx)) ' within bounds of current excerpt'];
  disp(logstr{end});
  logstr{end+1}=[int2str(length(find(ppIx))) ' positive peaks out of ' int2str(length(ppIx)) ' within bounds of current excerpt'];
  disp(logstr{end});

  npA=npA(npIx);
  npT=npT(npIx);
  ppT=ppT(ppIx);
  ppA=ppA(ppIx);

  % purge peaks which are too close to predecessor
  [npT,npIx]=tsldeadt(npT,minIPI);
  npA=npA(npIx);
  [ppT,ppIx]=tsldeadt(ppT,minIPI);
  ppA=ppA(ppIx);

  logstr{end+1}=[int2str(length(find(npIx))) ' negative peaks out of ' int2str(length(npIx)) ' with minimal IPI (' num2str(minIPI) ' ms)'];
  disp(logstr{end});
  logstr{end+1}=[int2str(length(find(ppIx))) ' positive peaks out of ' int2str(length(ppIx)) ' with minimal IPI (' num2str(minIPI) ' ms)'];
  disp(logstr{end});

  % now compute ipi
  npIPI=[diff(npT);nan];
  ppIPI=[diff(ppT);nan];

  % write local variables back to cell array..
  negPA{AP.LFPccInd(chaI)}=npA;
  negPT{AP.LFPccInd(chaI)}=npT;
  posPA{AP.LFPccInd(chaI)}=ppA;
  posPT{AP.LFPccInd(chaI)}=ppT;
  % AND fill preallocated cell arrays with IPI
  negPIPI{AP.LFPccInd(chaI)}=npIPI;
  posPIPI{AP.LFPccInd(chaI)}=ppIPI;

end

% generation of variables to be exported
if doExport
  for i=1:length(r)
    collPeak(i).npA=cell(AP.nAllLFPCh,1);
    collPeak(i).npIPI=cell(AP.nAllLFPCh,1);
    collPeak(i).ppA=cell(AP.nAllLFPCh,1);
    collPeak(i).ppIPI=cell(AP.nAllLFPCh,1);
  end
end
  
% --- second loop: loop over behaviors and compute CVs for each channel
for i=1:length(r)
  if ~isempty(r(i).iPts)
    disp([r(i).segmentType]);
    % preallocate arrays holding CVs in individual segments for all channels
    npCvA=repmat(nan,[r(i).ni  AP.nLFPCh]);
    npCvIPI=npCvA;
    ppCvA=npCvA;
    ppCvIPI=npCvA;

    % convert iPts to ms
    iMs=discrete2cont(r(i).iPts,WP.osi*.001,'intv',1);

    for chaI=1:AP.nLFPCh
      disp([rawCh(AP.LFPInd(chaI)).nm ' ..']);

      % assign to local variables again
      npA=negPA{AP.LFPccInd(chaI)};
      npT=negPT{AP.LFPccInd(chaI)};
      ppA=posPA{AP.LFPccInd(chaI)};
      ppT=posPT{AP.LFPccInd(chaI)};
      npIPI=negPIPI{AP.LFPccInd(chaI)};
      ppIPI=posPIPI{AP.LFPccInd(chaI)};

      % point to peaks inside interval defined by first and last segment
      % of current behavior
      npIx=npT>=iMs(1,1) & npT<=iMs(end,2);
      ppIx=ppT>=iMs(1,1) & ppT<=iMs(end,2);
      % omit last peak because ipi could be nan
      npIx(end)=0;
      ppIx(end)=0;

      npA=npA(npIx);
      npT=npT(npIx);
      npIPI=npIPI(npIx);
      ppA=ppA(ppIx);
      ppT=ppT(ppIx);
      ppIPI=ppIPI(ppIx);

      npIx=repmat(logical(1),size(npA));
      ppIx=repmat(logical(1),size(ppA));
      % linear index to first element not yet used
      npFix=1;
      ppFix=1;
      % compute CVs segment by segment
      for g=1:r(i).ni
        % --- negative peaks
        % update indices:
        tmpIx=npT(npIx)>=iMs(g,1);
        % note usage of find not compatible with V6
        npFix=npFix+find(tmpIx,1)-1;
        npIx(1:npFix-1)=0;
        % use var tmpIx to identify events < right border and to index npA,
        % npIPI and so on
        tmpIx=npIx;
        tmpIx(npFix:end)=tmpIx(npFix:end) & npT(npIx)<iMs(g,2);

        npCvA(g,chaI)=std(npA(tmpIx))/mean(npA(tmpIx));
        npCvIPI(g,chaI)=std(npIPI(tmpIx))/mean(npIPI(tmpIx));

        if doExport
          collPeak(i).npA{AP.LFPccInd(chaI)}=cat(1,collPeak(i).npA{AP.LFPccInd(chaI)},npA(tmpIx));
          collPeak(i).npIPI{AP.LFPccInd(chaI)}=cat(1,collPeak(i).npIPI{AP.LFPccInd(chaI)},npIPI(tmpIx));
        end

        % --- positive peaks
        % update indices:
        tmpIx=ppT(ppIx)>=iMs(g,1);
        % note usage of find not compatible with V6
        ppFix=ppFix+find(tmpIx,1)-1;
        ppIx(1:ppFix-1)=0;
        % use var tmpIx to identify events < right border and to index ppA,
        % ppIPI and so on
        tmpIx=ppIx;
        tmpIx(ppFix:end)=tmpIx(ppFix:end) & ppT(ppIx)<iMs(g,2);

        ppCvA(g,chaI)=std(ppA(tmpIx))/mean(ppA(tmpIx));
        ppCvIPI(g,chaI)=std(ppIPI(tmpIx))/mean(ppIPI(tmpIx));
        
        if doExport
          collPeak(i).ppA{AP.LFPccInd(chaI)}=cat(1,collPeak(i).ppA{AP.LFPccInd(chaI)},ppA(tmpIx));
          collPeak(i).ppIPI{AP.LFPccInd(chaI)}=cat(1,collPeak(i).ppIPI{AP.LFPccInd(chaI)},ppIPI(tmpIx));
        end
        
      end
    end % for:chaI=1:AP.nLFPCh
    % assign results to fields of r
    eval(['r(i).' STshort 'NegPeakCvA=npCvA;']);
    eval(['r(i).' STshort 'NegPeakCvAMn=mean(npCvA);']);
    eval(['r(i).' STshort 'NegPeakCvAStd=std(npCvA);']);
    eval(['r(i).' STshort 'NegPeakCvIPI=npCvIPI;']);
    eval(['r(i).' STshort 'NegPeakCvIPIMn=mean(npCvIPI);']);
    eval(['r(i).' STshort 'NegPeakCvIPIStd=std(npCvIPI);']);

    eval(['r(i).' STshort 'PosPeakCvA=ppCvA;']);
    eval(['r(i).' STshort 'PosPeakCvAMn=mean(ppCvA);']);
    eval(['r(i).' STshort 'PosPeakCvAStd=std(ppCvA);']);
    eval(['r(i).' STshort 'PosPeakCvIPI=ppCvIPI;']);
    eval(['r(i).' STshort 'PosPeakCvIPIMn=mean(ppCvIPI);']);
    eval(['r(i).' STshort 'PosPeakCvIPIStd=std(ppCvIPI);']);
  end % if:~isempty(r(i).iPts)
end % for:length(r)

% export?
if doExport
  save([AP.resPath '\' AP.resFn '_oscillReg'],'collPeak');
end