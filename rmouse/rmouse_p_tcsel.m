function rmouse_p_tcsel
% ** function rmouse_p_tcsel

global AP WP

% ------- plots of time course: selected parameters
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


fh6=mkfig('TimeCourse_sel'); 
orient landscape;
labelscale('fontSz',8,'scaleFac',1,'lineW',.75,'markSz',2.5); 
% figure title 
subplot('position',[WP.xmarg .95 1-2*WP.xmarg .02]);
set(gca,'xlim',[0 1],'ylim',[0 1]);
axis off
th=text(.5,.5,WP.figTitle,'fontsize',12,'fontweight','bold','horizontalalignment','center');
% vertical extent available to subplots
ve=.87;
% horizontal extent available to subplots
he=.94;
ph={[]};
sph={[]};

% determine contiguous segments so lines can be plotted
% elementary time step given excerpt length and overlap (the factor 
% accounts for rounding errors)
minTStep=discrete2cont(AP.ppSeg-AP.dftOlapPts,WP.osi*.001)/6e4*1.1;
for bi=1:length(ubix)
  tGroupIx{bi}=[0; find(diff(T(ubix{bi}))>minTStep); length(ubix{bi})];
end
  

% plot
for rvi=1:nrv
  sph{rvi}=subplot('position',[2*WP.xmarg  ve/nrv*(nrv-rvi)+WP.ymarg*1.8  he-2*WP.xmarg  ve/nrv-WP.ymarg*.15]);  
  hold on;
  for bi=1:length(ubix)
    pcol=AP.segmentType{ub(bi),3};
    psymb=AP.segmentType{ub(bi),4};    
%     ph=plot(T(ubix{bi}),Y(ubix{bi},rvi),psymb);
%     set(ph,'color',pcol);   
    for g=1:length(tGroupIx{bi})-1
      ix=tGroupIx{bi}(g)+1:tGroupIx{bi}(g+1);
      ph=plot(T(ubix{bi}(ix)),Y(ubix{bi}(ix),rvi),[psymb '-']);
      set(ph,'color',pcol);
    end
  end
  % parameters have to be treated separately
  [co,ncadh,bins]=cumh(Y(:,rvi),.001,'p',[.03 .97]);  
  subplot(sph{rvi}), 
  axis tight;
  xl=get(gca,'xlim');
  if isempty(co) | ~all(isfinite(co)) | diff(co)==0, 
    set(gca,'color',[.75 .75 .75],'xlim',xl);      
  else
    set(gca,'color',[.75 .75 .75],'ylim',co,'xlim',xl);
  end 
  grid on, box on
  ti=rv{rvi};  
  switch rv{rvi}
    case {'thgaeCCPeak','thgaeCCPeakT','thgaeCCPeakPhase',...
        'thgaeCCPosPeakDecay','thgaeCCZScore','thgaeCCZTestP',...
        'thPosPeakCvA','thPosPeakCvIPI','thNegPeakCvA','thNegPeakCvIPI'}
      ti=[ti ', within site'];      
    case {'rawDePE','rawThPE','rawThNarrowPE','rawBePE','rawGaPE',...
        'rawGaNarrowPEMn','rawRiPE','rawPPeak','rawPPeakT',...
        'thCCPosPeakDecay',...
        'gaePPeak','gaePPeakT',...
        'gaeThPEMn','gaeThNarrowPEMn'}
      ti=[ti ', auto'];            
    otherwise
      ti=[ti ', cross, {\Delta}=' num2str(dDist*1000) ' {\mu}m'];
  end
  ultext(ti,.05);
  if rvi<nrv
    set(gca,'xtick',[]);
  else
    xlabel('time (min)');
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
    print(pa,[WP.figName '_tc' ext]);   
  end
end