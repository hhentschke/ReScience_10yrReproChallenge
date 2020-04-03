% hypothesis: strong gamma bursts enhance theta amplitude and briefly reset (reduce) lags 
% between theta at s.p. and s.lm

% ideas: 
% - create phase plane for behavioral types separately & plot together
% - identify temporary 'phase resets'/abnormal phases on the basis of

colormap(jet);
nLev=30;
% IN 5 and IN 11
chInd1=6;
chInd2=12;
% excerpt
% mix
exc=[0 600];
% immobile
exc=[1442  1570];
% exploring
exc=[80 240];

excPts=cont2discrete(exc*1e6,osi,'intv',1);
% prealloc, taking into account downsampling fac
dfac=2;
d=repmat(0,floor((diff(excPts)+1)/dfac),2);
% load & downsample data:
% x axis: theta princ ch (putative s.l-m)
tmpd=strmread([AP.strmDir '\' rawCh(chInd2).thetaFn],'intv',excPts,'verbose',0);
d(:,1)=tmpd(1:dfac:end);
% y axis: theta (putative s.p.)
tmpd=strmread([AP.strmDir '\' rawCh(chInd1).thetaFn],'intv',excPts,'verbose',0);
d(:,2)=tmpd(1:dfac:end);

% dependent var: gamma envelope princ ch (putative s.l-m)
tmpd=strmread([AP.strmDir '\' rawCh(chInd1).gammaEnvFn],'intv',excPts,'verbose',0);            
dd(:,1)=tmpd(1:dfac:end);

% some test sines
% n=size(d,1);
% d(:,1)=((1:n)/n*100*pi)';
% d(:,2)=d(:,1)-pi*.8;
% d=sin(d)+.01*(rand(size(d))-.5);

% extend arrows' length by this factor
magFac=5;

% velocity (mV/s)
dv=diff(d,1,1)/(osi*1e-6);

[n1,n2]=size(d);
dmax=max(d,[],1);
dmin=min(d,[],1);
dspan=dmax-dmin;
nBins=50;
% bin borders for histc: nBins+1
binBs=repmat(((0:nBins)/nBins)',1,n2).*repmat(dspan,nBins+1,1)+repmat(dmin,nBins+1,1);
% find indices to events - dimension by dimension (binBs may be modified to contain
% unequal number of bin borders in the future)

% preallocate structure array of cell arrays holding indices to events in bins
dim(1).name='theta SLM';
dim(1).nBins=nBins;
dim(1).binBs=binBs(:,1);
% centers of bins
dim(1).binC=dim(1).binBs(1:end-1)+.5*diff(dim(1).binBs(1:2));
dim(1).ix=cell(1,dim(1).nBins);

dim(2).name='theta SP';
dim(2).nBins=nBins;
dim(2).binBs=binBs(:,2);
dim(2).binC=dim(2).binBs(1:end-1)+.5*diff(dim(2).binBs(1:2));
dim(2).ix=cell(1,dim(2).nBins);

arrTempl=repmat(nan,dim(1).nBins,dim(2).nBins);
% arrays holding averaged vector components in this dimension
dim(1).D=arrTempl
dim(2).D=arrTempl;
% array holding number of points in each coordinate
N=arrTempl;
% array holding some velocity (see below) in each coordinate
V=arrTempl;
% array holding gamma envelope magnitude in each coordinate
E=arrTempl;


% x and y arrays for plotting
[x,y]=ndgrid(dim(1).binC,dim(2).binC);

% get 1D indices
for i=1:length(dim)
  for k=1:nBins
    % restrict d because one parameter to be computed further down is vel
    dim(i).ix{k}=find(d(1:end-1,i)>=dim(i).binBs(k) & d(1:end-1,i)<dim(i).binBs(k+1));
  end
end

% for each n-tupel in parameter space compute the average n-dim vector
for i=1:dim(1).nBins
  for j=1:dim(2).nBins
    tmpix=intersect(dim(1).ix{i},dim(2).ix{j});
    N(i,j)=length(tmpix);
    if ~isempty(tmpix)
      if length(tmpix)==1
        dim(1).D(i,j)=dv(tmpix,1);
        dim(2).D(i,j)=dv(tmpix,2);
        % V(i,j)=sqrt(dv(tmpix,1).^2 + dv(tmpix,2).^2);        
        % velocity of SP theta
        V(i,j)=dv(tmpix,2);        
        E(i,j)=dd(tmpix);
      else
        dim(1).D(i,j)=mean(dv(tmpix,1));
        dim(2).D(i,j)=mean(dv(tmpix,2));
        % mean vel
        % V(i,j)=mean(sqrt(dv(tmpix,1).^2 + dv(tmpix,2).^2));
        V(i,j)=mean(dv(tmpix,2));
        E(i,j)=mean(dd(tmpix));        
      end
    else
      dim(1).D(i,j)=nan;
      dim(2).D(i,j)=nan;
      % zeroes look nicer on plots
      V(i,j)=nan;
      E(i,j)=0;      
    end
  end
end
% use sumfig..
clf

[c,cph]=contourf(x,y,E,nLev);
set(cph,'lineStyle','none');
nicexyax(10);
axlims=[get(gca,'xlim') get(gca,'ylim')];
hold on
% philosophy: first channel = first dimension in the matrices covering parameter
% space = row indices = x axis on plot. Consequently, both the 'D' matrices and the 
% y bins have to be flipped/rotated etc. in multiple ways
lh=quiver(dim(1).binC,flipud(dim(2).binC),flipud((dim(1).D)'),flipud((dim(2).D)'),magFac);
set(lh,'color','k');
axis(axlims);
xlabel(dim(1).name);
ylabel(dim(2).name);


break
orient landscape
subplot(2,2,1)
% philosophy: first channel = first dimension in the matrices covering parameter
% space = row indices = x axis on plot. Consequently, both the 'D' matrices and the 
% y bins have to be flipped/rotated etc. in multiple ways
lh=quiver(dim(1).binC,flipud(dim(2).binC),flipud((dim(1).D)'),flipud((dim(2).D)'),magFac);
xlabel(dim(1).name);
ylabel(dim(2).name);
nicexyax;
axlims=[get(gca,'xlim') get(gca,'ylim')];


% density of points: plot log with offset of 1
subplot(2,2,2);
[c,cph]=contourf(x,y,log(N+1),nLev);
set(cph,'lineStyle','none');
axis(axlims);
title('density of pts');

% vel
subplot(2,2,3);
[c,cph]=contourf(x,y,V,nLev);
set(cph,'lineStyle','none');
axis(axlims);
title('velocity');

% vel
subplot(2,2,4);
[c,cph]=contourf(x,y,E,nLev);
set(cph,'lineStyle','none');
axis(axlims);
title('av gamma env magn');

hold on
% axis manual
% for i=1:200, plot(d([1:20]+5*i,1),d([1:20]+5*i,2),'ro-'), pause(.05); drawnow; end
