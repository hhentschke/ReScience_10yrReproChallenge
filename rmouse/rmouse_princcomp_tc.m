function rmouse_princcomp_tc
% ** function rmouse_princcomp_tc
% computes principal components from combinations of segment-wise
% computed, selected parameters 

% improvements:
% - text display of variability

global AP WP

% choose parameters to collect - must be a field name of r (will be put in eval)
% the ones that seem to differ most between behaviors
rv={'rawPPeakT','rawPPeak',...
  'rawThPE','rawBePE','rawGaPE',...
  'thCCPeak','thCCPeakT','thNegPeakCvA',...
  'thgaeCCPeak','thgaeCCPeakT',...
};
% behaviors 
behav={'immobile','exploring'};
% identify channel 400 um more dorsal of principal (=lm) channel for CC plots
dDist=.4;

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


% how many PC should we take a look at?
nPC=4;
% what kinda plot?
plotmd='scatter';
plotmd='contour';

% delete columns with (i) ANY nans or (ii) zero variability
kix=find(any(isnan(Y),1) | ~any(diff(Y,1,1),1));
if ~isempty(kix)
  warning('the following results vars contain nans or have zero variance and cannot be used for PC:');
  disp(strvcat(rv{kix}));
end
Y(:,kix)=[];
% delete corresponding field names in rv
rv(kix)=[];
nYc=size(Y,2);
nrv=length(rv);
if nYc<nPC
  warning('cutting down number of principal components to number of parameters extracted');
  nPC=nYc;
end

% *** principal components
% normalize data (divide by std)
Y=Y./repmat(std(Y,1,1),size(Y,1),1);
[pcs,nd,variances,t2] = princomp(Y);
% pcs = the actual 'principal components'; each component sits in a ROW 
% nd = 'component scores', the data in the new coordinate system defined by the principal
% components (same size as original data)


% ************* plot *********************
fh7=mkfig('PrincComp_tc'); 
labelscale('scaleFac',1.0,'fontSz',8,'lineW',1.0,'markSz',2); 
orient landscape
colormap(flipud(gray));
% for each behavior plot 50 and 75% contour levels (normalized to max bin of each behavior) 
v=[.33 .66];
% number of combinations of PCs
[c,ncombin]=combin(nPC,'autoCC',0);
% depending on that, design layout of subplots (approx 3:4 ratio)
nRow=ceil(sqrt(3/4*ncombin));
nCol=ceil(ncombin/nRow);
% number of bin BORDERS per dim: max of .2*sqrt(number of points), 9
nBin=max(9,.2*sqrt(length(T)));

for ii=1:ncombin
  cmpi=c(ii,2);
  refi=c(ii,1);
  subplot(nRow,nCol,ii), hold on
  % first plot: contour or the like 
  % in setting axis limits ignore outlying 2%
  xl=(cumh(nd(:,refi),.01,'p',[.01 .99]))';
  yl=(cumh(nd(:,cmpi),.01,'p',[.01 .99]))';
  
  switch plotmd
    case 'contour'
      % n by n bins
      xb=linspace(xl(1),xl(2),nBin);
      yb=linspace(yl(1),yl(2),nBin);
      % the corresponding variants to be used as axis values for plots
      xb_pl=xb(1:end-1)+.5*diff(xb(1:2));
      yb_pl=yb(1:end-1)+.5*diff(yb(1:2));      
      h=hist2d([nd(:,cmpi) nd(:,refi)],yb,xb);
      ph=image(xb_pl,yb_pl,h);
      set(ph,'cdatamapping','scaled'); 
      axis([xl yl]);
      axis manual
      % overlay contours for individual behaviors
      for bi=1:length(ubix)
        pcol=AP.segmentType{ub(bi),3};      
        h=hist2d([nd(ubix{bi},cmpi) nd(ubix{bi},refi)],yb,xb);
        % normalize
        h=h/max(max(h));
        % contour requires a 'mock' line spec here because otherwise it creates patches
        % instead of lines 
        [cm,cph]=contour(xb_pl,yb_pl,h,v,'g');
        set(cph(:),'color',pcol,'linewidth',1.5);     
      end
    case 'scatter'
      for bi=1:length(ubix)      
        pcol=AP.segmentType{ub(bi),3};
        ph=plot(nd(ubix{bi},refi),nd(ubix{bi},cmpi),'s');
        set(ph,'color',pcol);
      end      
  end
  title(['PC ' int2str(refi) ' vs PC ' int2str(cmpi)]);
end
drawnow

% print?
if ~isempty(AP.printas{1}), 
  for i=1:length(AP.printas)
    pa=AP.printas{i};
    if strfind(pa,'ps'), ext='.ps';
    elseif strfind(pa,'jpeg'), ext='.jpg';
    else ext='';
    end
    print(pa,[WP.figName '_pc' ext]);   
  end
end