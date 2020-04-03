function abf2mat(job)
%  **function abf2mat(job)
% loads complete abf file, does sth to the data in it and saves the result
% as a mat file so that the data can be read by the rmouse matDload
% function. Needs DS as a global variable.
% The job may be either of the following:
% - 'resample' - resamples data (target sampling freq is 1000 Hz)
% - 'hifi' - hipass filters data at ~0.2 Hz, thus removing any base line 
% offset or very slow fluctuation

global DS WP

rmouse_ini;

% 'targeted sampling interval' to which data will be re-sampled 
tsi=1000;
ch=DS.rawCh(:,1);
stop=10;
stop='e';

% if dpath does not contain a drive letter, pre-pend WP.rootPath
if isempty(strfind(DS.dpath,':')), DS.dpath=[WP.rootPath DS.dpath]; end
savFn=[DS.dpath '\' DS.abfFn '.mat']

% delete any previously generated files
butt='yes';
if exist(savFn,'file')
  butt=questdlg('will overwrite previously generated file - proceed?');
end

if strcmpi(butt,'yes')
  for chInd=1:length(ch)
    [D,si]=abfload([DS.dpath '\' DS.abfFn '.abf'],'start',0,'stop',stop,'channels',ch(chInd));
    % ** disregard polarity (DS.rawSignalInverted) here because the matfile
    % produced here wil be treated in the same way as abf files in rmouse (i.e.
    % the raw data will be inverted, if necessary)
    % save channel by its deblanked name in DS
    dbch=ch{chInd};
    dbch=dbch(~(int16(dbch)==32));
    switch job
      case 'resample'
        D=resample(D,si,tsi);
        abfi.si=tsi;
      case 'hifi'
        D=hifi(D,si,0.2,'rs',20);
        abfi.si=si;
      otherwise
        error('bad job');
    end
    % convert to single
    D=single(D);
    % assign
    eval([dbch '=D;']);

    % save in matfile
    if chInd==1, 
      % put all relevant information needed in rmouse in structure abfi
      eval(['abfi.dataPtsPerChan=length(' dbch ');']);
      save(savFn,dbch,'abfi','-mat');
    else 
      save(savFn,dbch,'-mat','-append');
    end
    % delete
    eval(['clear ' dbch ';']);    
  end
end



