% creates summary plots of iso+atropine experiments 
% -> to run load results file produced by combine_r
rmouse_ini;
ornt='portrait';
figdir='c:\projects\rmouse\paper_atropine\rawFig\';
figName=mfilename;
printas=[]; '-dps2';
labelscale('fontSz',8,'scaleFac',.8,'lineW',0.8,'markSz',5);
figure(1), clf, orient(ornt);

rvix=strmatch('thNegPeakCvA',rv);
Rraw.d{1,rvix}=abs(Rraw.d{1,rvix});
Rraw.d{2,rvix}=abs(Rraw.d{2,rvix});


for spi=1:6;
  switch spi
    case 1
      curRv='rawPPeakTMn_auto' ;  curti='{\theta} peak frequency'; ylab='Hz';
      curRv='thNegPeakCvIPI';  curti='{\theta} CV_{IPI}'; ylab='';
    case 2
      curRv='rawPPeakMn_auto' ;  curti='{\theta} peak amplitude'; ylab='mV^2/Hz';
    case 3
      curRv='thNegPeakCvA';  curti='{\theta} CV_{Amplitude}'; ylab='';
    case 4
      curRv='rawGaPEMn_auto' ;  curti='{\gamma} power'; ylab='mV^2';
    case 5
      curRv='thgaeCCPeakMn_auto';  curti='C_{{\theta},{\gamma}Env}, peak correlation'; ylab='';
    case 6
      curRv='thgaeCCZScore_auto';  curti='C_{{\theta},{\gamma}Env}, Z score'; ylab='';
    otherwise
      return
  end

  subplot(4,3,spi), cla
  hold on

  rvix=strmatch(curRv,rv);
  % immobile
  icollP=[];
  pcIx=Rraw.d{1,rvix}(:,1)==0;
  icollP=(cat(1,icollP,Rraw.d{1,rvix}(pcIx,2:4)))';
  ph1=plot(1:3,icollP,'o');
  set(ph1,'color','k','markerfacecolor','k');
  tmpph=plot(2:3,icollP(2:3,:),'-');
  set(tmpph,'color','k');
  tmpph=plot(1:2,icollP(1:2,:),':');
  set(tmpph,'color',[.7 .7 .7]);
  
  % dirty ANOVA
  dd=icollP(:);
  dd=[dd repmat((1:3)',3,1) [1 1 1 2 2 2 3 3 3]'];
  rmaov1(dd,0.05)

  % small ttests on the side
  [h,p(1)]=ttest(icollP(1,:)',icollP(2,:)');
  [h,p(2)]=ttest(icollP(2,:)',icollP(3,:)');
  [h,p(3)]=ttest(icollP(1,:)',icollP(3,:)');
  [curti ': p(ctrl vs. iso)= ' num2str(p(1))...
    '; p(iso vs. iso+atr)= ' num2str(p(2))...
    '; p(ctrl vs. iso+atr)= ' num2str(p(3))]
  
  
  ecollP=[];
  pcIx=Rraw.d{2,rvix}(:,1)==0;
  ecollP=(cat(1,ecollP,Rraw.d{2,rvix}(pcIx,2:4)))';
  % plot exploring data only for control..
  ph2=plot(1,ecollP(1,:),'s');
  set(ph2,'color',[.5 .5 .5]);
  % but connect to 'immobile' data points for drug conditions
  tmpph=plot(1:2,[ecollP(1,:);icollP(2,:)],':');
  set(tmpph,'color',[.7 .7 .7]);

  set(gca,'xtick',1:3,'xticklabel',strvcat({'control','   iso   ','iso+atrop'}));
  nicexyax(10);
  % finally, put broken lines in background
  ch=get(gca,'children');
  bli=strmatch(':',get(ch,'LineStyle'),'exact');
  ch=[ch(setdiff(1:length(ch),bli)); ch(bli)];
  set(gca,'children',ch);
  
  title(curti);
  ylabel(ylab);
  if p<=.05, urtext('*',.95); end
    
  if spi==4, 
    legend([ph1(1) ph2(1)],'immobile','exploring');
  end

  % {'wt2141','wt2001','wt2206'}
end

if ~isempty(printas),
  print(printas,[figdir figName]);
end
