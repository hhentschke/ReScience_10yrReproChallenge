% --- analysis settings
curRv='rawgaeCohPeak_cross';

% --- plot settings
% print? 
printas='-dpsc2';
% printas=[];
% appearance of plots
labelscale('fontSz',10,'scaleFac',1,'lineW',.5,'markSz',10); 
% rmouse_ini;
ornt='landscape';
figdir='d:\projects\hippocampus\rmouse\paper_atropine\rawFig\';
figName=mfilename;
nColors=128;
figure(1), clf

% --- load data
load rcomb_wtAtropine_rawgaeCoh

% --- preliminaries
idv(1).name='rec site';
% rec sites; conventions: depth in mm; 0 is principal channel, spacing 0.1,
% values decreasing in dorsal direction
idv(1).level=[-0.6; -0.5; -0.4; -0.3; -0.2; -0.1; 0.0];
idv(1).nLevel=numel(idv(1).level);
% some shorties
behIx=strmatch('behavior',{RInfo.name});
druIx=strmatch('drug',{RInfo.name});
subIx=strmatch('subject',{RInfo.name});

nBehav=numel(RInfo(behIx).level);
nDrug=numel(RInfo(druIx).level);
nSubj=numel(RInfo(subIx).level);

rvIx=strmatch(curRv,rv,'exact');

% preallocate container for collected data
template_collD=repmat(nan,[idv(1).nLevel,nSubj,idv(1).nLevel]);

% --- extract & plot
for behi=1:nBehav
  for drui=1:nDrug
    tmpSubjIx=Rraw.indv{behi,rvIx,1}(:,1,1);
    collD=template_collD;
    for subi=1:nSubj
      % rec Sites
      tmp=Rraw.d{behi,rvIx,1}(tmpSubjIx==subi,1,1);
      [ux,ix,ox]=intersect(round(10*tmp),round(10*idv(1).level));
      % data proper
      tmp=Rraw.d{behi,rvIx,1}(tmpSubjIx==subi,drui+1,:);
      collD(ox,subi,ox)=tmp(ix,1,ix);
    end
    % --- plot
    subplot(nBehav,nDrug,(behi-1)*nDrug+drui);
    switch curRv
      case 'rawgaeCohPeak_cross'
        cLim=[0 0.75];
        colormap(flipud(coma('orangeblackblue','ncols',nColors)));
      case 'rawgaeCohPeakF_cross'
        cLim=[7 10];
        colormap(summer(nColors));
      case 'rawgaeCohTh_cross'
        cLim=[0 0.4];
        colormap(spring(nColors));
    end
    cohMat=permute(nanmean(collD,2),[1 3 2]);
    % cohMat=permute(nanstd(collD,0,2),[1 3 2]);
    ih=imagesc(cohMat,cLim);
    cph=gca;
    set(cph,'xaxisloc','top',...
      'ytick',[1:idv(1).nLevel],'xtick',[1:idv(1).nLevel]);
    axis square
    xlabel('gammaEnv'); % sic!
    ylabel('raw');
    cbh=colorbar;
%     if ci==n2 && ri==n1
%       cbh=colorbar;
%     end

    
  end
end
      



if ~isempty(printas),
  print(printas,[figdir figName]);
end

