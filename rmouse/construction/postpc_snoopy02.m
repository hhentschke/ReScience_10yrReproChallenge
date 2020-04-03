function helper02(intv)
% ** function helper02(intv)
% an open script probing this and that:
% - compute CC between two channels from one large stretch of data, starting with one
% (=the whole) segment, proceeding to do this for an increasing number of segments, and
% compute the average deviation of the individual segments' CC from the large one


global DS AP WP 

% minimal and maximal length of subinterval (s)
minLenSubInt=0.2;
maxLenSubInt=diff(intv);
% cc lag
ccLag=.2
% RELATIVE overlap of segments
rOlapSubInt=.4;
% length of subintervals (s)
base=2;
lenSubInt=minLenSubInt*base.^[log2((maxLenSubInt-.002)/minLenSubInt):-.5:0]';
nDiv=size(lenSubInt,1);
% where to fish for peaks (s)
peakLag=.035;

ccScaleOpt='coeff';'coeff_ub';

rmouse_ini;
rmouse_APcheck;
rawCh=rmouse_chan;

strmType={'delta','theta','thetaEnv','thetaLo','thetaLoEnv','thetaHi','thetaHiEnv','gamma','gammaEnv'};
strmType={'theta'};
nStrms=length(strmType);

if isempty(strfind(DS.dpath,':')), DS.dpath=[WP.rootPath DS.dpath]; end
if isempty(strfind(AP.strmDir,':')), AP.strmDir=[WP.rootPath AP.strmDir]; end

% get si
abfi=abfinfo([DS.dpath '\' DS.abfFn '.abf'],'verbose',0);      
si=abfi.si;

% --- convert all relevant vars to points
ccLag_pts=cont2discrete(ccLag*1e6,si,'intv',1);
intv_pts=cont2discrete(intv*1e6,si,'intv',1);
tP=cont2discrete(peakLag,si*1e-6,'intv',0)+ccLag_pts;

  
% generate one var for each stream, this is most flexible for all sorts of plots,
% including pllplot overlays
for ci=length(AP.LFPccInd):-1:1
  for six=1:length(strmType)
    % sic: rawCh(ci)
    eval(['d(:,ci)=strmread([AP.strmDir ''\'' rawCh(ci).' strmType{six} 'Fn],''intv'',intv_pts,''verbose'',0);' ]);
  end
end


% t=(1:size(d,1))'/1000*2*pi;
% 
% d(:,1)=sin(t*8);
% d(:,2)=sin(t*8+.52);


close all
figure
hold on
% cc for whole segment
CC0=zeros(2*ccLag_pts+1,1);
ph1=plot(CC0);
% mean deviation
ccMnDev=zeros(nDiv,1);
for g=1:nDiv
  [subInt,subInt_pts,nilen,nilen_pts]=mkintrvls([0 diff(intv)],'resol',si*1e-6,'ilen',lenSubInt(g),...
    'olap',lenSubInt(g)*rOlapSubInt,'border','skip');
  nSubInt=size(subInt_pts,1);
  ccDev=zeros(1,nSubInt);
  cc=zeros(2*ccLag_pts+1,nSubInt);
  delete(ph1);
  for subi=1:nSubInt
    % - cc  
    [cc(:,subi),lag]=xxcorr(d(subInt_pts(subi,1):subInt_pts(subi,2),1),...
      d(subInt_pts(subi,1):subInt_pts(subi,2),2),...
      ccLag_pts,ccScaleOpt);
%     if g==1 & subi==1
%       CC0=cc;
%     else
%       ccDev(subi)=sqrt(mean((CC0-cc(:,subi)).^2));
%     end
  end
  CC0=mean(cc,2);
  %  ccDev=sqrt(mean((repmat(CC0,1,nSubInt)-cc).^2));
  tmpr=evdeal(cc,'idx',{'closepeak'},'tP',tP);
  ph1=plot(lag,cc,'k-');
  ph0=plot(lag,CC0,'r-');
  set(ph0,'linewidth',2.5);
  drawnow;
  % ccMnDev(g)=mean(ccDev);
  ccMnDev(g)=std(tmpr.tlPosPeakT);
end

figure
subplot(3,1,1)
ph=plot(lenSubInt,ccMnDev,'o-');
% set(gca,'xscale','log','yscale','log');
axis tight
subplot(3,1,2)
ph=plot(lenSubInt,ccMnDev,'o-');
set(gca,'xscale','log');
axis tight
subplot(3,1,3)
ph=plot(lenSubInt,ccMnDev,'o-');
set(gca,'xscale','log','yscale','log');
axis tight
xlabel('win length (s)');
ylabel('mean distance');
    










% counterpart to strmwrite
function d=strmread(varargin)
global DS
d=i16load(varargin{:})/2^15*DS.nsRng;
