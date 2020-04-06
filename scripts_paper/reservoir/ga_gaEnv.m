% Generates plot of gamma power vs theta power of gammaEnv in atropine 
% condition. Both parameters are normalized to their control value. The
% idea is to estimate to what degree the theta-modulated component of gamma 
% is reduced (compared to total gamma) by atropine
% --- graphics
labelscale('fontSz',8,'scaleFac',.275,'lineW',.6,'markSz',4); 
ornt='portrait';
printas='-dps2';
printas=[];
figDir='d:\projects\hippocampus\rmouse\paper_atropine\rawFig\';
figName=mfilename;

load rcomb_wtAtropine


% parameters to be plotted against each other
curRv{1}='rawGaPEMn_auto' ;  axl{1}='gamma power';
curRv{2}='gaeThPEMn_auto' ;  axl{2}='theta power of gamma Env';

% curRv{1}='gaeThPEMn_auto' ;  axl{1}='theta power of gamma Env';
% curRv{2}='thgaeCCPeakMn_auto' ;  axl{1}='thgaeCCPeakMn';


nIndv=5;
% individuals down the columns, column order:
% 1st par ctrl | 2nd par ctrl | 1st par atr | 2nd par atr
% behaviors in slices
%  (note: for plots below we need only the normalized atropine data;
%  however, the structure of collP as outlined above is maintained for
%  possible extensions of this plot in the future)
collP=repmat(nan,[nIndv,4,2]);
% loop over parameters
for g=1:2
  % find parameter
  rvix=strmatch(curRv{g},rv);
  % loop over behaviors
  for bi=1:2
    % identify princ chan
    pcIx=Rraw.d{bi,rvix}(:,1)==0;
    collP(:,[g g+2],bi)=Rraw.d{bi,rvix}(pcIx,[2 3]);
  end
end

% % normalize collP to control values
% collP(:,[1 3],:)=collP(:,[1 3],:)./repmat(collP(:,1,:),[1 2 1]);
% collP(:,[2 4],:)=collP(:,[2 4],:)./repmat(collP(:,2,:),[1 2 1]);

% multiply gammaEnv power by 1000 so axis labels look better
collP(:,[2 4],:)=collP(:,[2 4],:)*1000;

figure(1), clf, hold on
% Jan 2007: skip data for immobility because the figure will be overcrowded
if 0
  % --- immobility
  bi=1;
  ph0=plot(collP(:,[1 3],bi)',collP(:,[2 4],bi)','ko-');
  ph1=plot(collP(:,1,bi),collP(:,2,bi),'ko');
  ph2=plot(collP(:,3,bi),collP(:,4,bi),'ko');
  set(ph2,'markerfacecolor','k')
  % % regression
  % idepD=[ones(nIndv,1) collP(:,3,bi)];
  % alpha=.05;
  % [rcoeff,confI,tmpnix,tmpnix2,stats]=regress(collP(:,4,bi),idepD,alpha);
  % disp(['y=' num2str(rcoeff(2),2) '*x+' num2str(rcoeff(1),2)...
  %   '; r^2=' num2str(stats(1)) '; p=' num2str(stats(3))]);
end
% --- exploring
bi=2;
ph3=plot(collP(:,[1 3],bi)',collP(:,[2 4],bi)','ks:');
ph4=plot(collP(:,1,bi),collP(:,2,bi),'ks');
ph5=plot(collP(:,3,bi),collP(:,4,bi),'ks');
set(ph5,'markerfacecolor','k')
% % regression 
% idepD=[ones(nIndv,1) collP(:,3,bi)];
% alpha=.05;
% [rcoeff,confI,tmpnix,tmpnix2,stats]=regress(collP(:,4,bi),idepD,alpha);
% disp(['y=' num2str(rcoeff(2),2) '*x+' num2str(rcoeff(1),2)...
%   '; r^2=' num2str(stats(1)) '; p=' num2str(stats(3))]);        

% axis([.2 1 .2 1])
nicexyax;
xlabel(curRv{1});
ylabel(curRv{2});

% line(get(gca,'xlim'),get(gca,'ylim'),'linestyle','--','color','k');
if exist('ph0','var')
  legend([ph1 ph2 ph4 ph5],...
    {'im, ctrl','im, atr','ex, ctrl','ex, atr'},'location','SouthEast')
else
  legend([ph4 ph5],...
    {'ctrl','atropine'},'location','SouthEast')
end
% % finally, paired ttesting of similarity of reduction of both parameters
% we need to normalize collP to control values
collP(:,[1 3],:)=collP(:,[1 3],:)./repmat(collP(:,1,:),[1 2 1]);
collP(:,[2 4],:)=collP(:,[2 4],:)./repmat(collP(:,2,:),[1 2 1]);
bi=1;
[h,p]=ttest(collP(:,3,bi),collP(:,4,bi))
bi=2;
[h,p]=ttest(collP(:,3,bi),collP(:,4,bi))


if ~isempty(printas), print(printas,[figDir figName]); end

