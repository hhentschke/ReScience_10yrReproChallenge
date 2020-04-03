function abfmerge2mat_rmouse(varargin)
% ** function abfmerge2mat_rmouse(varargin) merges abf data files, heeding
% timing, and saves data to matfile. This is a variant working with data
% structures used by rmouse.m, so global DSET and ANPAR are expected.
% !! It is assumed that neither channels nor sampling interval nor any 
% other viable variable changes from one recording to the next !!
%                         >>> INPUT VARIABLES >>>
%
% NAME        TYPE/DEFAULT             DESCRIPTION
% ch          cell arr or 'all'/'all'  channels to be read & merged
% savFn       char arr, <name of 1st   name of merged mat file
%             file in list> 
%                         <<< OUTPUT VARIABLES <<<
%
% NAME        TYPE/DEFAULT             DESCRIPTION
%
%


global DSET ANPAR
% DSET and ANPAR must be 1D and have at least 2 elements in chronological order
nF=length(DSET);

rmouse_ini;

% --- defaults
ch=DSET(1).rawCh(:,1);
savFn=[WP.rootPath DSET(1).dpath '\' DSET(1).abfFn '.mat'];
pvpmod(varargin)
if ischar(ch) & strcmpi(ch,'all')
  % all's fine (hopefully)
end

stop='e';

% determine length of merged streams
[nix,nix2,abfi1]=abfload([DSET(1).dpath '\' DSET(1).abfFn '.abf'],'info');
[nix,nix2,abfi2]=abfload([DSET(end).dpath '\' DSET(end).abfFn '.abf'],'info');
% abfi1=abfinfo([DSET(1).dpath '\' DSET(1).abfFn '.abf']);
% abfi2=abfinfo([DSET(end).dpath '\' DSET(end).abfFn '.abf']);

% abfi is a 'made-up' structure containing only vital information about the
% merged file needed by rmouse
abfi.si=abfi1.si;
% .recTime is specified in seconds
abfi.dataPtsPerChan=cont2discrete(abfi2.recTime(2)-abfi1.recTime(1),abfi.si*1e-6,'intv',1);
% this timing information is needed by combine_bv
abfi.lFileStartTime=abfi1.lFileStartTime;
% preallocate
D=repmat(0,abfi.dataPtsPerChan,1);

% delete any previously generated files
butt='yes';
if exist(savFn,'file')
  butt=questdlg('will overwrite previously generated file - proceed?');
end


if strcmpi(butt,'yes')
  % do things channel by channel, then file by file (less memory demand than file by file)
  for chInd=1:length(ch)
    dbch=ch{chInd};
    disp(['channel: ' dbch '...']);
    % deblank names
    dbch=dbch(~(int16(dbch)==32));
    for fi=nF:-1:1
      [nix,nix2,abfiCurr]=abfload([DSET(fi).dpath '\' DSET(fi).abfFn '.abf'],'info');
      % abfiCurr=abfinfo([DSET(fi).dpath '\' DSET(fi).abfFn '.abf']);
      % allocate right slots..
      ix=cont2discrete(abfiCurr.recTime-abfi1.recTime(1),abfi.si*1e-6,'intv',1);
      % ** load: disregard polarity (DS.rawSignalInverted) here because the matfile
      % produced here wil be treated in the same way as abf files in rmouse (i.e.
      % the raw data will be inverted, if necessary)
      [tmpD,si]=abfload([DSET(fi).dpath '\' DSET(fi).abfFn '.abf'],'start',0,'stop',stop,'channels',ch(chInd));
      switch diff(ix)+1-length(tmpD)
        case 0
          ;
        case -1
          disp('real number of data points one more than expected');
          tmpD(end)=[];
        otherwise
          error('sampling time index not OK');
      end

      % embed
      D(ix(1):ix(2))=tmpD;
    end
    % save channel by its deblanked name in DS
    eval([dbch '=single(D);']);
    % save in matfile
    if chInd==1, 
      save(savFn,dbch,'abfi','-mat');
    else 
      save(savFn,dbch,'-mat','-append');
    end
    % delete
    eval(['clear ' dbch ';']);    
  end
end



