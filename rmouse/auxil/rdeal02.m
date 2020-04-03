% This script does a session-by-session comparison of selected results
% variables as computed by rmouse.m and combined for all sessions by 
% combine4b.m. Comparison is done for control vs. drug; it consists of
% curve fitting +running an F-test (like in combine4b.m). To run this
% script, load the contents of auto.mat or cross.mat. Either of these mat 
% files contains a results variable R.
% R is the collection of results from several experiments, won by running
% combine4b.m. It is a struct with (currently) 6 fields. Each field is a n-D cell
% array:
% - row=behavior (order as in behav)
% - col=parameter (order as in rv)
% - slice=genotype (1st Wt, 2nd mutant)
% Description of the fields (i.e., each individual cell's content):
%         d: collected data: 2d arr, 1st col electrode pos, 2nd col value (drug
%            exprmnts: 3rd+ cols = values ar var. concentrations)
%      indv: code for individual animal/session
%      ueix: 1d cell array; for each of the ue, these are the indices into the
%            corresponding R.d
%        ga: 'grand average': 2d arr, holding  mean|std|N
%     bstat:
%        ue:



% mm, limits of electrode depth (inclusive; slm=0, dorsal ones negative, ventral ones positive)
depthLim=[-.601 0];
% assume
% - all parameters were analyzed for each animal/conc/behav
uIndv=unique(R.indv{1,1});
% colors and symbols for plotting
pset={'b','m';'o','s'};
% multiple small figs
labelscale('fontSz',8,'scaleFac',1,'lineW',2,'markSz',5);

% switch q
%   case 'auto'
%     % all that make sense as autos 
%     rv_commented={...
%       'detheCCPeakMn', 'delta & theta envelope, peak crosscorrelation, dimensionless';...
%       'detheCCPeakTMn', 'delta & theta envelope, lag at peak crosscorrelation, ms';...
%       'thgaeCCPeakMn', 'theta & gamma envelope, peak crosscorrelation, dimensionless';...
%       'thgaeCCPeakTMn', 'theta & gamma envelope, lag at peak crosscorrelation, ms';...
%       'thgaeCCZScore', 'theta & gamma envelope, Z-score of peak crosscorrelation, dimensionless';...
%       'thgaeCCZTestP', 'theta & gamma envelope, p-value of peak crosscorrelation, dimensionless';...
%       'rawPPeakMn', 'theta, peaks of individual (segmental) power spectra, averaged, mV^2/Hz';...
%       'rawPPeakTMn', 'theta, frequency of peaks of individual (segmental) power spectra, Hz';...
%       'rawPMnPeak', 'theta, peak of averaged power spectra, mV^2/Hz';...
%       'rawPMnPeakT', 'theta, frequency peak of averaged power spectra, mV^2/Hz';...
%       'rawDePEMn', 'delta, power, mV^2';...
%       'rawThPEMn', 'theta, power, mV^2';...
%       'rawThNarrowPEMn', 'narrow band of theta (peak freq +/-1 Hz), power, mV^2';...
%       'rawBePEMn', 'beta, power, mV^2';...
%       'rawGaPEMn', 'gamma, power, mV^2';...
%       'rawRiPEMn', 'ripples, power, mV^2';...
%       'thCCPosPeakDecayMn', 'theta, height of largest side peak in autocorrelation, dimensionless';...
%       'thgaeCCPosPeakDecayMn', 'theta & gamma envelope, height of largest side peak in their correlation, dimensionless';...
%       'gaePPeakMn', 'gamma envelope, peaks of individual (segmental) power spectra, averaged, mV^2/Hz';...
%       'gaePPeakTMn', 'gamma envelope, frequency of peaks of individual (segmental) power spectra, Hz';...
%       'gaePMnPeak', 'gamma envelope, peak of averaged power spectra, mV^2/Hz';...
%       'gaePMnPeakT', 'gamma envelope, frequency peak of averaged power spectra, mV^2/Hz';...
%       'gaeThPEMn', 'gamma envelope, power in theta band, mV^2';...
%       'gaeThNarrowPEMn', 'gamma envelope, power in narrow theta band, mV^2';...
%       };
%     % currently of interest
%     rvix='all';
%     %rvix=[9];
% 
%     
%   case 'cross'
%     % all that make sense as cross-measures
%     rv_commented={...
%       'deCCPeakMn', 'delta, peak crosscorrelation, dimensionless';...
%       'deCCPeakTMn', 'delta, lag of peak crosscorrelation, ms';...
%       'thCCPeakMn', 'theta, peak crosscorrelation, dimensionless';...
%       'thCCPeakTMn', 'theta, lag of peak crosscorrelation, ms';...
%       'gaCCPeakMn', 'gamma, peak crosscorrelation, dimensionless';...
%       'gaCCPeakTMn', 'gamma, lag of peak crosscorrelation, ms';...
%       'thHieCCPeakMn', 'theta hi, peak crosscorrelation, dimensionless';...
%       'thHieCCPeakTMn', 'theta hi, lag of peak crosscorrelation, ms';...
%       'thLoeCCPeakMn', 'theta lo, peak crosscorrelation, dimensionless';...
%       'thLoeCCPeakTMn', 'theta lo, lag of peak crosscorrelation, ms';...
%       'gaeCCPeakMn', 'gamma envelope, peak crosscorrelation, dimensionless';...
%       'gaeCCPeakTMn', 'gamma envelope, lag of peak crosscorrelation, ms';...
%       'rawPPeakMn', 'theta, peaks of individual (segmental) power spectra, averaged, mV^2/Hz';...
%       'rawPPeakTMn', 'theta, frequency of peaks of individual (segmental) power spectra, Hz';...
%       'rawPMnPeak', 'theta, peak of averaged power spectra, mV^2/Hz';...
%       'rawPMnPeakT', 'theta, frequency peak of averaged power spectra, mV^2/Hz';...
%       'rawCohMnDe', 'delta, coherence, dimensionless';...
%       'rawCohMnTh', 'theta, coherence, dimensionless';...
%       'rawCohMnThNarrow', 'narrow band of theta, coherence, dimensionless';...
%       'rawCohMnBe', 'beta, coherence, dimensionless';...
%       'rawCohMnGa', 'gamma, coherence, dimensionless';...
%       'rawCohMnRi', 'ripples, coherence, dimensionless';...
%       'rawDePEMn', 'delta, power, mV^2';...
%       'rawThPEMn', 'theta, power, mV^2';...
%       'rawThNarrowPEMn', 'narrow band of theta (peak freq +/-1 Hz), power, mV^2';...
%       'rawBePEMn', 'beta, power, mV^2';...
%       'rawGaPEMn', 'gamma, power, mV^2';...
%       'rawRiPEMn', 'ripples, power, mV^2';...
%       'gaeCohMnTh', 'gamma envelope, coherence in theta range, dimensionless';...
%       'gaeCohMnThNarrow', 'gamma envelope, coherence in narrow theta range, dimensionless';...
%       'gaeCohMnDe', 'gamma envelope, coherence in delta range, dimensionless';...
%       };


for rvInd=1:length(rv)
  proceed=0;
  if strcmpi(q,'cross')
    % if ~isempty(strmatch(rv{rvInd},{'thCCPeakMn','thCCPeakTMn','gaCCPeakMn','gaCCPeakTMn','rawCohMnDe','rawCohMnTh','rawCohMnTh_narrow','rawCohMnBe','rawCohMnGa','gaeCohMnTh','gaeCohMnTh_narrow','gaeCohMnDe'},'exact'))
    if ~isempty(strmatch(rv{rvInd},{'thCCPeakTMn'},'exact'))
    % if ~isempty(strmatch(rv{rvInd},{'thCCPeakPhaseMn'},'exact'))
      proceed=1;
    end
  else
    % if ~isempty(strmatch(rv{rvInd},{'thgaeCCPeakTMn'},'exact'))
    % if ~isempty(strmatch(rv{rvInd},{'rawPMnPeak','rawThPEMn','rawGaPEMn','gaePMnPeak'},'exact'))      
    % if ~isempty(strmatch(rv{rvInd},{'rawGaNarrowPEMn'},'exact'))
    % if ~isempty(strmatch(rv{rvInd},{'thCCPosPeakDecayMn'},'exact'))
    % if ~isempty(strmatch(rv{rvInd},{'thgaeCCZTestP'},'exact'))
    if ~isempty(strmatch(rv{rvInd},{'thgaeCCZScore'},'exact'))    
    % if ~isempty(strmatch(rv{rvInd},{'thgaeCCPeakTMn'},'exact'))
    % if ~isempty(strmatch(rv{rvInd},{'thgaeCCPeakPhaseMn'},'exact'))
    % if ~isempty(strmatch(rv{rvInd},{'rawGaPEMn'},'exact'))
    % if ~isempty(strmatch(rv{rvInd},{'thCCPeakDecayMn'},'exact'))
      proceed=1;
    end
  end

  if proceed
    figure(rvInd), clf,
    set(gcf,'pos',[10   143   723   813]);
    for iInd=1:length(uIndv)
      try
        for bInd=1:length(behav)
          % identify data belonging to current session
          csInd=find(R.indv{bInd,rvInd}==uIndv(iInd));
          % ds1 ctrl, ds2 drug
          ds1=R.d{bInd,rvInd}(csInd,[1 2]);
          ds2=R.d{bInd,rvInd}(csInd,[1 3]);
          % restrict depth range
          ds1(ds1(:,1)<depthLim(1) | ds1(:,1)>depthLim(2),:)=[];
          ds2(ds2(:,1)<depthLim(1) | ds2(:,1)>depthLim(2),:)=[];
          % invert sign
          ds1(:,1)=ds1(:,1)*-1;
          ds2(:,1)=ds2(:,1)*-1;
          % sort
          ds1=sortrows(ds1,1);
          ds2=sortrows(ds2,1);
          % - combination of both
          ds12=sortrows([ds1; ds2]);

          % depending on parameter transform data and set up model
          [ft_,fo_,st_,ds1ix,ds2ix,ds12ix]=curveFit2rmousePar(ds1,ds2,ds12,rv{rvInd});

          % fit and determine quality of fit for wt and ko
          set(fo_,'Startpoint',st_);
          ds1f=fit(ds1(ds1ix,1),ds1(ds1ix,2),ft_ ,fo_);
          ds2f=fit(ds2(ds2ix,1),ds2(ds2ix,2),ft_ ,fo_);
          ds12f=fit(ds12(ds12ix,1),ds12(ds12ix,2),ft_ ,fo_);

          % --- part 4: create curves representing the fits & plot
          % 1. all data pts & fit
          fitx=ds12(1,1):(ds12(end,1)-ds12(1,1))/200:ds12(end,1);
          % cfit objects will be 'fevaluated' automatically
          ds1fit=ds1f(fitx);
          ds2fit=ds2f(fitx);
          ds12fit=ds12f(fitx);

          % plot
          titl=[behav{bInd} ', ' rv{rvInd} ', ' q];
          fnp=['WTdrugE_' behav{bInd} '_' rv{rvInd} '_' q];

          subplot(length(uIndv),length(behav),2*(iInd-1)+bInd), hold on
          title(titl);
          % first
          ph=plot(ds1(:,1),ds1(:,2),pset{2,1});
          set(ph,'color',pset{1,1});
          ph=plot(fitx,ds1fit,'-');
          set(ph,'color',pset{1,1});
          % second
          ph=plot(ds2(:,1),ds2(:,2),pset{2,2});
          set(ph,'color',pset{1,2});
          ph=plot(fitx,ds2fit,'-');
          set(ph,'color',pset{1,2});
          % combo fit in black
          plot(fitx,ds12fit,'k-');
          nicexyax;
          drawnow

          % --- part 5: F-test
          % statistical test for similarity (order: wt crf, ko crf, combo crf)
          [p,F,radj1,radj2]=curvecomp([ds1(ds1ix,:) ds1f(ds1(ds1ix,1))], ...
            [ds2(ds2ix,:) ds2f(ds2(ds2ix,1))],...
            [ds12(ds12ix,:) ds12f(ds12(ds12ix,1))],...
            length(ds1ix)-length(st_),length(ds2ix)-length(st_));

          urtext(['p=' num2str(p,'%1.3f')],.85,'fontsize',12);
          if p<.05
            disp(['***** H0 (identity of depth-response profiles) rejected, p= ' num2str(p)]);
            %       if p<.01, urtext('**',.9,'fontsize',20);
            %       else urtext('*',.9,'fontsize',20);
            %       end
          else
            disp(['H0 (identity of depth-response profiles) not rejected, p= ' num2str(p)]);
            %       urtext('n.s.');
          end
          % information on goodness of fit
          ultext(['R_{adj}(1):' num2str(radj1,3) '; R_{adj}(2):' num2str(radj2,3)]);
          ds1f
          ds2f
          ds12f
        end
      catch
        disp('current data set produced an error')
      end
    end
    pause
  end % if:proceed
end