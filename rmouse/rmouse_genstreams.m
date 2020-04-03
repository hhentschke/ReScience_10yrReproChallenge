function rmouse_genstreams
% ** function rmouse_genstreams
% generates 'streams': take one channel at a time, filter it, do the
% Hilbert transforms, and save all streams as int16 to individual binary
% files
% ** ALL data of one channel are digested, irrespective of the settings of
% AP.rawExcerpt (unless quickNdirty~=0)

global DS AP WP

quickNdirty=0;
if quickNdirty
  tmpstop=WP.eoe*.001;
else
  tmpstop='e';
end

if ~exist(AP.strmDir,'dir'),
  pdi=strfind(AP.strmDir,'\');
  mkdir(AP.strmDir(1:pdi(end)-1),AP.strmDir(pdi(end)+1:end));
end

% ----- loop over channels
for i=1:AP.nCh
  if exist([DS.dpath '\' DS.abfFn '.mat'],'file')
    [d,si]=matDload([DS.dpath '\' DS.abfFn '.mat'],'start',0,'stop',tmpstop,'channels',{AP.rawChAnNm{i}});
  elseif exist([DS.dpath '\' DS.abfFn '.raw'],'file')
    rawFInfo=rawinfo('filename',[DS.dpath '\' DS.abfFn '.raw'],'print','no');
    abfi=rawfi2genfi(rawFInfo,'channels',AP.rawChAnNm(i));
    si=abfi.si;
    if ischar(tmpstop) & strcmp(tmpstop,'e')
      d=rawload([DS.dpath '\' DS.abfFn '.raw'],AP.rawChAnNm(i),'full',rawFInfo);
    else
      d=rawload([DS.dpath '\' DS.abfFn '.raw'],AP.rawChAnNm(i),[0 tmpstop*1000],rawFInfo);
    end
    % put into array and convert to mV
    d=cat(2,d{:})/1000;
  else
    [d,si]=abfload([DS.dpath '\' DS.abfFn '.abf'],'start',0,'stop',tmpstop,'channels',{AP.rawChAnNm{i}});
  end
  [n1,n2]=size(d);
  if DS.rawSignalInverted
    d=-1*d;
  end

  strmList=AP.strm;
  nStrm=numel(strmList);
  sct=0;
  % ----- loop over requested streams
  % (do this using as little memory as possible, that is, for each channel,
  % one stream after the other, deleting the new streams as soon as they're
  % saved to disk)
  while ~isempty(strmList) || sct<nStrm
    sct=sct+1;
    fiFreq=eval(['AP.' strmList{1} 'CFreq']);
    strmFn=[DS.abfFn '_' AP.rawChAnNm{i} '_' strmList{1} '.i16'];
    filtD=bafi(d,si,fiFreq,'rs',AP.rs);
    % get rid of line hum?
    if ismember(strmList{1},{'gammaNarrow','gamma'}) && ~isempty(AP.lineFreq)
      filtD=killhum(filtD,si,AP.lineFreq);
    end
    % save
    try
      fid=fopen([AP.strmDir '\' strmFn],'w','ieee-le');
    catch
      % if it did not work, probably the directory does not exist, so write files
      % in current dir
      disp(['writing stream failed: ' lasterr '; now saving in current directory']);
      fid=fopen(strmFn,'w','ieee-le');
    end
    c=strmwrite(fid,filtD);
    fclose(fid);
    if c~=n1, 
      error(['sth went wrong writing to ' strmFn]); 
    end
  
    % envelope?
    if numel(strmList)>1 && strcmp([strmList{1} 'Env'],strmList{2})
      sct=sct+1;
      strmFn=[DS.abfFn '_' AP.rawChAnNm{i} '_' strmList{2} '.i16'];
      filtDEnv=abs(hilbert(filtD));
      % save
      try
        fid=fopen([AP.strmDir '\' strmFn],'w','ieee-le');
      catch
        disp(['writing stream failed: ' lasterr '; now saving in current directory']);
        fid=fopen(strmFn,'w','ieee-le');
      end
      c=strmwrite(fid,filtDEnv);
      fclose(fid);
      % clear stream
      filtDEnv=[];
      if c~=n1,
        error(['sth went wrong writing to ' strmFn]);
      end
      % remove names of written streams from list
      strmList(1:2)=[];
    else
      % remove name of written stream from list
      strmList(1)=[];
    end
    
  end
  
end % for: i=1:AP.nCh
