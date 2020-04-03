function combine_snoopy01
% ** function combine_snoopy01
% - combines individual rmouse data sets and plots time course of selected variables
% - is a blend of rmouse_p_tcSel and combine2 
% - needs ANPAR and DSET (=concatenation of individual and matching AP and DS) as global vars:
%   first row=control, second+ rows = post-drug-administration files, 
%   with only one column! For computation of grand averages or the like see combine3, combine3b 
% - the only (optional) input argument is an array of points in time (unit: minutes) 
%   where markers should be placed on the uppermost subplot

global ANPAR DSET

[n1 n2 n3]=size(ANPAR);
[dn1 dn2 dn3]=size(DSET);
if ~isequal([n1 n2 n3],[dn1 dn2 dn3]), error('DSET and ANPAR must have equal dimension'); end
if n1>1 | n3>1, error('DSET and ANPAR must have no more than 1 row'); end

% behaviors that should make it on the plots
behav={'immobile','exploring'};
nBehav=length(behav);
% corresponding colors and plot symbols
pCol={[.6 .4 .1],[.5 1 .5]};
pSymb={'s','o'};
% choose parameters to collect and plot - must be a field name of r
% (will be put in eval)
rv={'thCCPeak','thCCPeakT'};
rv={'thgaeCCPeak','thgaeCCPeakT'}
rv={'gaeCCPeak','gaeCCPeakT'}
nrv=length(rv);
% corresponding bins for histograms
binArr={(-1:.01:1)', (-50:1:100)'};

% cell array of histograms: each cell one parameter; 
% each element: one hist per column, as many columns as files, behaviors in
% slices
Yhist=cell(1,n2);
for ii=1:nrv
  Yhist{ii}=zeros(length(binArr{ii})-1,n2,nBehav);
end

% identify channel 500 um more dorsal of principal (=lm) channel for CC plots
dDist=.6;
% variable collecting parameters: first nrv cols=data, 
% two columns after that: global time | behavType
Y=[];
% unique numerical codes for all behaviors actually collected
ub=[];

rmouse_ini;

printas={[]};{'-djpeg90'};
figName=[WP.rootPath ANPAR(1).resPath '\' DSET(1).abfFn '_' DSET(end).abfFn];

% -------- PART I: collection of data
% loop over data sets for one experiment
for ii=1:n2
  AP=ANPAR(ii);
  DS=DSET(ii);
  % if paths do not contain a drive letter, pre-pend WP.rootPath
  if isempty(strfind(DS.dpath,':')), DS.dpath=[WP.rootPath DS.dpath]; end
  if isempty(strfind(AP.resPath,':')), AP.resPath=[WP.rootPath AP.resPath]; end    
  % extract si and timing information from abf file - if matfile exists, pick it instead of abf file
  if exist([DS.dpath '\' DS.abfFn '.mat'],'file')
    load([DS.dpath '\' DS.abfFn '.mat'],'abfi');      
  elseif exist([DS.dpath '\' DS.abfFn '.abf'],'file')
    abfi=abfinfo([DS.dpath '\' DS.abfFn '.abf'],'verbose',0);      
    % put out start and stop times in seconds from midnight
    disp(['recording start||end: ' sprintf('%6.3f || %6.3f',abfi.recTime) ' s from midnight']);
  else
    error([DS.dpath '\' DS.abfFn ' does not exist'])
  end
  % original sampling interval
  osi=abfi.si;
  % load results var..
  load([AP.resPath '\' AP.resFn],'r');
  % channel business
  rawCh=rmouse_chan;
  [doff,ccRefChInd]=min(abs(WP.elx(AP.LFPccInd)+dDist));
  ccRefChIdx=AP.LFPIdx(ccRefChInd);
  dDist=diff(WP.elx([ccRefChIdx AP.pcIdx]));
  % index to last row of current Y
  lix=size(Y,1);
  % time offset to add to local time
  if lix>0
    tOffs=abfi.lFileStartTime/60+discrete2cont(round(AP.ppSeg/2),osi*.001)/6e4;
  else 
    tOffs=abfi.lFileStartTime/60;
  end
  if ii==1
    figTitle=DS.aName;
  end
  % loop over behaviors, collect parameters
  for ri=1:length(r)
    lbi=strmatch(r(ri).segmentType(:,1),behav);
    if ~isempty(r(ri).iPts) & ~isempty(lbi)
      % local time axis - minutes, please
      t=discrete2cont(round(mean(r(ri).iPts,2)),osi*.001)/6e4;
      nt=length(t);
      Y=[Y; repmat(nan,nt,nrv+2)];
      % global time
      Y(end-nt+1:end,nrv+1)=t+tOffs;
      % code for behavior = index in local var behav
      Y(end-nt+1:end,nrv+2)=lbi;
      ub=[ub lbi];
      for rvi=1:nrv
        eval(['y=r(ri).' rv{rvi} ';']);
        % select extraction method depending on results var chosen 
        switch rv{rvi}
          case {'thgaeCCPeak','thgaeCCPeakT'}
            if ~isempty(y)
              y=y(:,AP.LFPpcInd1);
            else
              y=[];
            end
          case {'rawThPE','rawGaPE','rawPPeakMn','rawPPeakTMn'}
            if length(diag(y))>1
              y=y{AP.LFPpcInd2,AP.LFPpcInd2};
            else
              y=[];
            end
          otherwise
            if length(diag(y))>1
              y=y{ccRefChIdx,AP.LFPpcInd2};
            else
              y=[];
            end
        end % switch
        if isempty(y)
          warning(['data to be extracted do not exist']);
          % current rows of Y are preallocated with nans, so y can be left empty
        else
          % collect values
          Y(end-nt+1:end,rvi)=makecol(y);
          tmph=histc(y,binArr{rvi});
          Yhist{rvi}(:,ii,lbi)=tmph(1:end-1);
        end
      end % for:rv
    end % if:~isempty r.ipts & behav requested
  end % for:behaviors
end
% sort time stamps 
Y=sortrows(Y,nrv+1);
% subtract time offset
Y(:,nrv+1)=Y(:,nrv+1)-Y(1,nrv+1);
ub=unique(ub);



% ---------- PART II: plot -------------------------------------------------

for rvi=1:nrv
  % for each var one fig
  fh=mkfig(rv{rvi});
  orient tall;
  labelscale('fontSz',12,'scaleFac',1,'lineW',.25,'markSz',2);
  ct=0;
  for ii=1:n2
    for bi=1:nBehav
      ct=ct+1;
      subplot(n2,nBehav,ct);
      bar(binArr{rvi}(1:end-1),Yhist{rvi}(:,ii,bi),1.0,'k')
    end
  end
  % expand plots a little
  % rexy('xfac',1,'yfac',1.1);
end


if ~isempty(printas{1}),
  for i=1:length(printas)
    pa=printas{i};
    if strfind(pa,'ps'), ext='.ps';
    elseif strfind(pa,'jpeg'), ext='.jpg';
    else ext='';
    end
    print(pa,[figName '_tc' ext]);
  end
end

% -------- local func ----------
function figha=mkfig(ftag)
figha=findobj('tag',ftag);
if isempty(figha), figha=figure;
else  figure(figha);
end
tmpScrSz=get(0,'Screensize');
tmpScrSz=tmpScrSz*.7+.1*rand;
tmpScrSz(1)=tmpScrSz(1)+diff(tmpScrSz([1 3]))*.10+.01*rand;
tmpScrSz(2)=tmpScrSz(2)+diff(tmpScrSz([2 4]))*.10+.01*rand;  
set(figha,'position',tmpScrSz,...
  'tag',ftag,...
  'name',ftag,...
  'color',[.9 .9 1],...
  'numbertitle','off');
clf;
