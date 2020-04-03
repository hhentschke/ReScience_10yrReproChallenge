% ------- computes correlations of combinations of segment-wise computed, 
% selected parameters 
logstr{end+1}='---------- computing & plotting correlation coefficients of select variables..';

% extract parameters
rmouse_segParExtract;

% delete columns with ANY nans
kix=any(isnan(Y),1);
if ~isempty(kix)
  warning('the following results vars contain nans and cannot be used for PC:');
  disp(strvcat(rv{kix}));
end
Y(:,kix)=[];
% delete corresponding field names in rv
rv(kix)=[];
nYc=size(Y,2);
nrv=length(rv);


% ***** correlate (behaviors separately) ******
Rr=[];
Pp=[];
Ss=[];
for bi=1:length(ubix)
  [rr,pp]=corrcoef(Y(ubix{bi},:));
  Rr=cat(3,Rr,rr);
  Pp=cat(3,Pp,pp);
  % locate all non-redundant ones with significant correlation
  Ss=cat(3,Ss,triu(pp<=.05));
end
% collapse Ss, looking for significant combinations in either behavior
Ss=any(Ss,3);
% exclude trivial/uninteresting combinations
% Ss([1:2 5:12],:)=0;
% Ss(:,[1:2 5:12])=0;

[scrix,sccix]=find(Ss);
% how many?
nSC=length(scrix);

% ************* plot *********************
fh8=mkfig('corr_tc'); 
labelscale('scaleFac',1.0,'fontSz',8,'lineW',1.0,'markSz',2); 
orient landscape

% design layout of subplots (approx 3:4 ratio)
nRow=ceil(3/4*sqrt(nSC));
nCol=ceil(nSC/nRow);

for ii=1:nSC
  cmpi=scrix(ii);
  refi=sccix(ii);
  subplot(nRow,nCol,ii), hold on
  % in setting axis limits ignore outlying 2%
  xl=(cumh(Y(:,refi),.01,'p',[.01 .99]))';
  yl=(cumh(Y(:,cmpi),.01,'p',[.01 .99]))';

  for bi=1:length(ubix)
    pcol=AP.segmentType{ub(bi),3};
    ph=plot(Y(ubix{bi},refi),Y(ubix{bi},cmpi),'s');
    set(ph,'color',pcol);
  end
  axis([xl yl]);
  text(xl(1)+diff(xl)*.1,yl(2)-diff(yl)*.1,...
  strvcat(['r=' sprintf('%1.2f|%1.2f',squeeze(Rr(refi,cmpi,:)))],...
  ['p=' sprintf('%1.3f|%1.3f',squeeze(Pp(refi,cmpi,:)))]),...
  'fontweight','bold');
  title([rv{cmpi} ' vs ' rv{refi}]);
end

% print?
if ~isempty(AP.printas{1}), 
  for i=1:length(AP.printas)
    pa=AP.printas{i};
    if strfind(pa,'ps'), ext='.ps';
    elseif strfind(pa,'jpeg'), ext='.jpg';
    else ext='';
    end
    print(pa,[figName '_tcCorr' ext]);   
  end
end

clear Y T B pcs nd variances t2 ub ubix vex xb* yb* bi ii cmpi refi cm cph nBin nCol
