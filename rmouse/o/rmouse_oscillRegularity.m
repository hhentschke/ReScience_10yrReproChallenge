% This routine computes coefficients of variation (cv) of peak amplitude
% and inter-peak-intervals (ipi) of an oscillatory input stream (so far
% theta only) The idea is to find out how regular the oscillation is

% minimally acceptable ipi
minIPI=50;

% graphics
fh_oscReg=mkfig('oscillRegul');
orient tall
labelscale('fontSz',7,'scaleFac',1.0,'lineW',.5,'markSz',4);

% arrays for preallocation
% CVs 
nanny=repmat(nan,nLFPCh,1);
% histogram bins (100 each)
nBin=100;
binA=(-1.485:0.03:1.5)';
binIPI=(minIPI:(250-minIPI)/(nBin-1):250)';
% histograms
nanny2=repmat(nan,100,nLFPCh);
iCollect=[];

% load theta peaks (time is ms)
load(WP.thetaPFn,'negPA','negPT','posPA','posPT');

% delete this soon
for i=1:length(r)
  if ~isempty(r(i).iPts)
    r(i).firstPt=r(i).iPts(1,1);
  end
end

for chaI=1:nLFPCh
  % disp(['CV theta peak amp and ipi, ch ' rawCh(AP.LFPInd(chaI)).nm ' ..']);
  % assign to local variables
  npA=negPA{AP.LFPccInd(chaI)};
  npT=negPT{AP.LFPccInd(chaI)};
  ppA=posPA{AP.LFPccInd(chaI)};
  ppT=posPT{AP.LFPccInd(chaI)};

  % kick out all peaks which are not within bounds defined by AP.rawExcerpt 
  % - this speeds up computation of values for small excerpts
  mima=discrete2cont([min([r(:).firstPt]) max([r(:).lastPt])],osi*.001);
  
  npIx=npT>=mima(1) & npT<mima(2);
  ppIx=ppT>=mima(1) & ppT<mima(2);

  npA=npA(npIx);
  npT=npT(npIx);
  ppT=ppT(ppIx);
  ppA=ppA(ppIx);
  
  % purge peaks which are too close to predecessor
  [npT,npIx]=tsldeadt(npT,minIPI);
  npA=npA(npIx);
  [ppT,ppIx]=tsldeadt(ppT,minIPI);
  ppA=ppA(ppIx);

  % now compute ipi
  npIPI=[diff(npT);nan];
  ppIPI=[diff(ppT);nan];
  
  for i=1:length(r)
    if ~isempty(r(i).iPts)
      % disp([r(i).segmentType '..']);
      % preallocate
      if chaI==1
        r(i).thPosPeakCvA=nanny;
        r(i).thPosPeakCvIPI=nanny;
        r(i).thNegPeakCvA=nanny;
        r(i).thNegPeakCvIPI=nanny;
        iCollect=[iCollect i];
        tmpr(i).thNegPeakAH=nanny2;
        tmpr(i).thNegPeakIPIH=nanny2;
        tmpr(i).thPosPeakAH=nanny2;
        tmpr(i).thPosPeakIPIH=nanny2;
      end
      % convert ePts to ms
      eMs=discrete2cont(r(i).ePts,osi*.001,'intv',1);

      npIx=zeros(size(npT));      
      ppIx=zeros(size(ppT));
      
      % collect indices to peaks excerpt by excerpt  
      for g=1:r(i).ne
        npIx=npIx | (npT>=eMs(g,1) & npT<eMs(g,2));
        ppIx=ppIx | (ppT>=eMs(g,1) & ppT<eMs(g,2));
      end

      % finally, do the calculations
      r(i).thNegPeakCvA(chaI)=std(npA(npIx))/mean(npA(npIx));
      % for ipi, set last index to 0 because the last ipi is nan
      npIx(end)=0;
      r(i).thNegPeakCvIPI(chaI)=std(npIPI(npIx))/mean(npIPI(npIx));
      r(i).thPosPeakCvA(chaI)=std(ppA(ppIx))/mean(ppA(ppIx));
      ppIx(end)=0;
      r(i).thPosPeakCvIPI(chaI)=std(ppIPI(ppIx))/mean(ppIPI(ppIx));
      
      % normalized histogram data for plots
      nNp=sum(npIx);
      nPp=sum(ppIx);
      tmpr(i).thNegPeakAH(1:nBin,chaI)=hist(npA(npIx),binA)/nNp;
      tmpr(i).thNegPeakIPIH(1:nBin,chaI)=hist(npIPI(npIx),binIPI)/nNp;
      tmpr(i).thPosPeakAH(1:nBin,chaI)=hist(ppA(ppIx),binA)/nPp;
      tmpr(i).thPosPeakIPIH(1:nBin,chaI)=hist(ppIPI(ppIx),binIPI)/nPp;
      

      % export raw data for princ chan?
      if chaI==AP.LFPpcInd1
        enpA=npA(npIx);
        eppA=ppA(ppIx);
        enpIPI=npIPI(npIx);
        eppIPI=ppIPI(ppIx);
        save([AP.resPath '\' AP.resFn '_oscillReg_prinChan' r(i).segmentType],'enpA','eppA','enpIPI','eppIPI');
        clear enpA eppA enpIPI eppIPI
      end
      
    end % if:~isempty(r(i).iPts)
  end % for:length(r)
end

% export?
% save([AP.resPath '\' AP.resFn '_oscillReg'],'tmpr','binA','binIPI');

% plot - amplitudes left column, ipi right column
spc=0;
% restrict to expl and imm
iCollect=intersect([immV explV],iCollect);
tmpNegPeakAH=-1*cat(3,tmpr(iCollect).thNegPeakAH);
tmpNegPeakIPIH=-1*cat(3,tmpr(iCollect).thNegPeakIPIH);
tmpPosPeakAH=cat(3,tmpr(iCollect).thPosPeakAH);
tmpPosPeakIPIH=cat(3,tmpr(iCollect).thPosPeakIPIH);

for chaI=1:nLFPCh
  spc=spc+1;
  subplot(nLFPCh,3,spc); hold on;
  ph=stairs(binA,permute(tmpNegPeakAH(:,chaI,:),[1 3 2]));
  set(ph,{'color'},AP.segmentType(iCollect,3));
  ph=stairs(binA,permute(tmpPosPeakAH(:,chaI,:),[1 3 2]));
  set(ph,{'color'},AP.segmentType(iCollect,3));
  niceyax
  % xlabel('amplitude')
  ultext(DS.rawCh{AP.LFPIdx(chaI),1},.1);
 
  spc=spc+1;
  sph(spc)=subplot(nLFPCh,3,spc); hold on;
  ph=stairs(binIPI,permute(tmpNegPeakIPIH(:,chaI,:),[1 3 2]));
  set(ph,{'color'},AP.segmentType(iCollect,3));
  ph=stairs(binIPI,permute(tmpPosPeakIPIH(:,chaI,:),[1 3 2]));
  set(ph,{'color'},AP.segmentType(iCollect,3));
  niceyax
  % xlabel('IPI')
  ultext(DS.rawCh{AP.LFPIdx(chaI),1},.1);
  
  spc=spc+1;

end
drawnow

if ~isempty(AP.printas{1}),
  for i=1:length(AP.printas)
    pa=AP.printas{i};
    if strfind(pa,'ps'), ext='.ps';
    elseif strfind(pa,'jpeg'), ext='.jpg';
    else ext='';
    end
    print(pa,[figName '_thetaOscillReg' ext]);
  end
end

if ishandle(fh_oscReg), delete(fh_oscReg); end

% reset graphics to default
labelscale('fontSz',7,'scaleFac',1.0,'lineW',1.0,'markSz',4);

clear minIPI fh_oscReg nBin nann* tmp* binA binIPI iCollect spc sph npA npT ppA ppT npIx ppIx npIPI ppIPI nNp pNp