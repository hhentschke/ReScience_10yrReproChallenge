function rmouse_parcorr_tc
% ** function rmouse_parcorr_tc
% computes correlation between segment-wise values of different results
% parameters, the aim being to identify parameters that covary

global AP WP

% ------- selected parameters
% choose parameters to collect - must be a field name of r (will be put in eval)
rv={'thgaeCCPeak',...
  'rawGaPE',...  'gamma power'
  'rawThNarrowPE',...
  'gaeThNarrowPE',...
  'thCCPeak',...
  'gaeCCPeak',...
  'thNegPeakCvIPI',...
  'thNegPeakCvA',...
  'gaePosPeakCvIPI',...
  'gaePosPeakCvA',...
  };


% behaviors 
behav={'immobile','exploring'};
behav={'exploring'};
% identify channel 200 um more dorsal of principal (=lm) channel for CC plots
dDist=.5;

% extract parameters
[T,Y,B,ub,ubix]=rmouse_segparextract(rv,behav,dDist);

% delete columns of ALL nans and adjust parameters accordingly 
kix=all(isnan(Y),1);
if ~isempty(kix)
  Y(:,kix)=[];
  % delete corresponding field names in rv
  rv(kix)=[];
end
% number of results variables to be plotted 
nrv=numel(rv);


% compute corrs
m=corr(Y);

% plot #1: color coded corr matrix
fh6=mkfig('ParCorr_sel_m'); 
colormap(coma('bluered'));
orient landscape;
labelscale('fontSz',8,'scaleFac',1,'lineW',.75,'markSz',2.5); 
imagesc(m,[-1 1]);
colorbar
set(gca,'XAxisLocation','top');
set(gca,'ytick',1:nrv,'yticklabel',rv,'xtick',1:nrv);

% plot #2: scatter plots
fh7=mkfig('ParCorr_sel_scatter'); 
orient landscape;
[h,axH]=gplotmatrix(Y,[],B);
for g=1:nrv
  for h=1:nrv
    subplot(axH(h,g)),
    smarttext(rv{h},.02,.9);
    smarttext(rv{g},.02,.1);
  end
end

% print?
if ~isempty(AP.printas{1}), 
  for i=1:length(AP.printas)
    pa=AP.printas{i};
    if strfind(pa,'ps'), ext='.ps';
    elseif strfind(pa,'jpeg'), ext='.jpg';
    else ext='';
    end
    figure(fh6)
    print(pa,[WP.figName '_parCorr_m' ext]);   
    figure(fh7)
    print(pa,[WP.figName '_parCorr_scatter' ext]);   
  end
end

% restore standard fig settings
labelscale('fontSz',7,'scaleFac',1.0,'lineW',1.0,'markSz',4);