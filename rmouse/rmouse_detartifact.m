function af_etsl=rmouse_detartifact(artifV,af_tsl_ext)

global DS AP WP

etslconst;
af_etsl=[];

if exist([DS.dpath '\' DS.abfFn '.mat'],'file')
  [d,si]=matDload([DS.dpath '\' DS.abfFn '.mat'],'start',WP.boe*.001,'stop',WP.eoe*.001,'channels',AP.rawChMonNm);
elseif exist([DS.dpath '\' DS.abfFn '.raw'],'file')
  rawFInfo=rawinfo('filename',[DS.dpath '\' DS.abfFn '.raw'],'print','no');
  abfi=rawfi2genfi(rawFInfo,'channels',AP.rawChMonNm);  
  [d,si]=rawload([DS.dpath '\' DS.abfFn '.raw'],AP.rawChMonNm,[WP.boe WP.eoe],rawFInfo);
  d=cat(2,d{:})/1000;
  abfi=rawfi2genfi(rawFInfo,'channels',AP.rawChMonNm);
  si=abfi.si;
else
  [d,si]=abfload([DS.dpath '\' DS.abfFn '.abf'],'start',WP.boe*.001,'stop',WP.eoe*.001,'channels',AP.rawChMonNm);
end
[n1,n2]=size(d);
if DS.rawSignalInverted
  d=-1*d;
end
% lopass filter?
if ~isempty(AP.afCFreq),
  d=lofi(d,si,AP.afCFreq);
end
% threshold
tmptsl=tcd(d,si,AP.afThresh);
% convert to arrays and combine
af_tsl=[];
for chInd=1:n2
  af_tsl=cat(1,af_tsl,tmptsl{chInd});
end

if isempty(af_tsl) && isempty(af_tsl_ext) 
  disp('no artifacts detected');
  figure(WP.sumFigH);
  subplot('position',[WP.xmarg WP.ymarg .33/2-2*WP.xmarg .8/3*2-2*WP.ymarg]);
  text(.1,.5,'no artifacts');
  drawnow
else
  if ~isempty(af_tsl)
    % set time frame: t0=start of neural recording, see above
    af_tsl=af_tsl+WP.boe;
  end
  % *** if external artifact tsl is nonempty add it here, after the
  % internally detected artifacts have been time-shifted
  if ~isempty(af_tsl_ext)
    af_tsl=cat(1,af_tsl,af_tsl_ext);
  end

  % generate list of forbidden intervals:
  % 1. replace each event by two new events, namely at AP.afWin(1) and
  % AP.afWin(2)
  af_tsl=cat(1,af_tsl+AP.afWin(1)*1000,af_tsl+AP.afWin(2)*1000);
  % 2. sort
  af_tsl=sort(af_tsl);
  % 3. crop events beyond borders
  af_tsl(af_tsl-WP.boe<0 | WP.eoe-af_tsl<0)=[];
  % 4. get finished etsl from etslburstf
  % dead time in ms - one si must be added
  tmp=diff(AP.afWin)*1e3+si/1000;
  af_etsl=etslburstf(af_tsl,tmp,'recLen',WP.eoe);
  % 5. set correct tag=segment type
  af_etsl(:,etslc.tagCol)=artifV;
  msg=[int2str(length(af_tsl)) ' artifact events; their combination resulted in '...
    int2str(size(af_etsl,1)) ' bad periods'];
  disp(msg);
  % plot up to 20 data segments of length 2*AP.afWin around artifacts
  % (divide by number of channels)
  figure(WP.sumFigH);
  subplot('position',[WP.xmarg WP.ymarg .33/2-2*WP.xmarg .8/3*2-2*WP.ymarg]);
  naf=size(af_etsl,1);
  nseg=20;
  % for raw data excerpt d t0 is WP.boe>0, not WP.bor=0, so subtract
  % WP.boe
  ptsl=af_etsl(1:max([1,floor(length(AP.rawChMonNm)*naf/nseg)]):naf,etslc.tsCol)-WP.boe;
  pPts=cont2discrete([ptsl ptsl]+repmat(AP.afWin*2*1e3,size(ptsl,1),1),si*.001);
  % crop to existing number of points
  pPts(pPts(:,1)<1,1)=1;
  pPts(pPts(:,2)>n1,2)=n1;
  % collect in cell array (traces may have different lengths)
  for i=size(pPts,1):-1:1
    for chInd=length(AP.rawChMonNm):-1:1
      pD{(i-1)*length(AP.rawChMonNm)+chInd}=d(pPts(i,1):pPts(i,2),chInd);
    end
  end
  pllplot(pD,'si',si,'spacing','fixed','dy',5);
  th=title('artifacts'); set(th,'fontsize',11,'fontweight','bold');
  drawnow
end
