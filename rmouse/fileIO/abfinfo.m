function h=abfinfo(fn,varargin)
% ** function h=abfinfo(fn,varargin)
% puts out information about data files in the pclamp abf-format.
% Optional input parameters listed below (= all except the file name) 
% must be specified as parameter/value pairs, e.g. as in 
%          s=abfinfo('d:\data01.abf','verbose',0);
%
%                    >>> INPUT VARIABLES >>>
%
% NAME            TYPE/DEFAULT    DESCRIPTION
% fn              string          abf data file 
% machineF        char array,     the 'machineformat' input parameter of the
%                  'ieee-le'       matlab fopen function. 'ieee-le' is the correct 
%                                  option for windows; depending on the
%                                  platform the data were recorded/shall be read
%                                  by abfload 'ieee-be' is the alternative.
% verbose         scalar,1        if nonzero, all information will be 
%                                  printed on screen
%
%                    <<< OUTPUT VARIABLES <<<
%
% NAME            TYPE/DEFAULT    DESCRIPTION
% h            structure       contains most vital header parameters 
%                                  +some derived parameters with self-explanatory
%                                  field names

machineF='ieee-le';
verbose=1;
pvpmod(varargin);

if verbose
  disp(['** ' mfilename]);
end

% vital header parameters - initialize with -1
% name, position in header in bytes, type, value)
tmp=repmat(-1,1,16);
headPar={
  'fFileVersionNumber',4,'float',-1;  
  'nOperationMode',8,'int16',-1; 
  'lActualAcqLength',10,'int32',-1;
  'nNumPointsIgnored',14,'int16',-1;
  'lActualEpisodes',16,'int32',-1;
  'lFileStartTime',24,'int32',-1;
  'lDataSectionPtr',40,'int32',-1;
  'lSynchArrayPtr',92,'int32',-1;
  'lSynchArraySize',96,'int32',-1; 
  'nDataFormat',100,'int16',-1;            
  'nADCNumChannels', 120, 'int16', -1;
  'fADCSampleInterval',122,'float', -1; 
  'fSynchTimeUnit',130,'float',-1;
  'lNumSamplesPerEpisode',138,'int32',-1;        
  'lPreTriggerSamples',142,'int32',-1;        
  'lEpisodesPerRun',146,'int32',-1;        
  'fADCRange', 244, 'float', -1;
  'lADCResolution', 252, 'int32', -1;
  'nFileStartMillisecs', 366, 'int16', -1;
  'nADCPtoLChannelMap', 378, 'int16', tmp;
  'nADCSamplingSeq', 410, 'int16',  tmp;
  'sADCChannelName',442, 'uchar', repmat(tmp,1,10);
  'fADCProgrammableGain', 730, 'float', tmp;
  'fInstrumentScaleFactor', 922, 'float', tmp;
  'fInstrumentOffset', 986, 'float', tmp;
  'fSignalGain', 1050, 'float', tmp;
  'fSignalOffset', 1114, 'float', tmp;
  'nTelegraphEnable',4512,'int16',tmp;
  'fTelegraphAdditGain',4576,'float',tmp
};

fields={'name','offs','numType','value'};
s=cell2struct(headPar,fields,2);
numOfParams=size(s,1);

if ~exist(fn,'file'), error(['could not find file ' fn]); end
if verbose, disp(['opening ' fn '..']); end
[fid,messg]=fopen(fn,'r',machineF);
if fid == -1,error(messg);end;   
% determine absolute file size
fseek(fid,0,'eof');
fileSz=ftell(fid);
fseek(fid,0,'bof');

% read all vital information in header
% read value from header and put in structure h
for g=1:numOfParams,
  if fseek(fid, s(g).offs,'bof')~=0, 
    fclose(fid);
    error(['something went wrong locating ' s(g).name]); 
  end;
  sz=length(s(g).value);
  eval(['[h.' s(g).name ',n]=fread(fid,sz,''' s(g).numType ''');']);
  if n~=sz, error(['something went wrong reading value(s) for ' s(g).name]); end;
end;

if h.lActualAcqLength<h.nADCNumChannels, warning('less data points than sampled channels in file'); end;
% the numerical value of all recorded channels (numbers 0..15)
sc=h.nADCSamplingSeq(1:h.nADCNumChannels);

% tell me where the data start
BLOCKSIZE=512;
switch h.nDataFormat
case 0
  dataSz=2;  %bytes/point
  precision='int16';
case 1
  dataSz=4;  %bytes/point
  precision='float32';
otherwise
  error('invalid number format');
end;
h.headOffset=h.lDataSectionPtr*BLOCKSIZE+h.nNumPointsIgnored*dataSz;
% fADCSampleInterval is the TOTAL sampling interval
h.si=h.fADCSampleInterval*h.nADCNumChannels;

% channel information:
% the numerical value of all recorded channels (numbers 0..15)
% ! line below does not work with some abf files in which the nADCSamplingSeq array has values 0 instead of -1 for non-recorded channels!
% recChIdx=nADCSamplingSeq(nADCSamplingSeq~=-1);
recChIdx=h.nADCSamplingSeq(1:h.nADCNumChannels);
% the names of ALL 16 channels, e.g. 'IN 8'
recChNames=[reshape(char(h.sADCChannelName),10,16)]';
% recorded channels
for i=1:length(recChIdx)
  h.recChNames{i}=deblank(recChNames(recChIdx(i)+1,:));
end
if verbose
  disp('available channels:');
  disp(h.recChNames);
  disp(['sampling interval (us): ' num2str(h.si)]);  
end

switch h.nOperationMode
  case 1
    if verbose, disp('data were acquired in event-driven variable-length mode');   end
    % h.nSweeps=lActualEpisodes; % check whether that is correct
    % see abfload for parameters to be included into h
  case {2,5}
    % see abfload for further parameters to be included into h
    if verbose
      if h.nOperationMode==2, disp('data were acquired in event-driven fixed-length mode');
      else disp('data were acquired in waveform fixed-length mode (clampex only)');
      end;
    end

    % extract timing information on sweeps
    if (h.lSynchArrayPtr<=0 || h.lSynchArraySize<=0),
      fclose(fid);
      error('internal variables ''lSynchArraynnn'' are zero or negative');
    end
    switch h.fSynchTimeUnit
      case 0  % time information in synch array section is in terms of ticks
        h.synchArrTimeBase=1;
      otherwise % time information in synch array section is in terms of usec
        h.synchArrTimeBase=h.fSynchTimeUnit;
    end;

    % the byte offset at which the SynchArraySection starts
    h.lSynchArrayPtrByte=BLOCKSIZE*h.lSynchArrayPtr;
    % before reading Synch Arr parameters check if file is big enough to hold them
    % 4 bytes/long, 2 values per episode (start and length)
    if h.lSynchArrayPtrByte+2*4*h.lSynchArraySize>fileSz,
      fclose(fid);
      error('file seems not to contain complete Synch Array Section');
    end
    
    if fseek(fid,h.lSynchArrayPtrByte,'bof')~=0,
      fclose(fid);
      error('something went wrong positioning file pointer to Synch Array Section');
    end
    
    [synchArr,n]=fread(fid,h.lSynchArraySize*2,'int32');
    if n~=h.lSynchArraySize*2,
      fclose(fid);
      error('something went wrong reading synch array section');
    end
    
    % make synchArr a h.lSynchArraySize x 2 matrix
    synchArr=permute(reshape(synchArr',2,h.lSynchArraySize),[2 1]);
    if numel(unique(synchArr(:,2)))>1
      fclose(fid);
      error('sweeps of unequal length in file recorded in fixed-length mode');
    end
    % the length of sweeps in sample points (**note: parameter lLength of
    % the ABF synch section is expressed in samples (ticks) whereas
    % parameter lStart is given in synchArrTimeBase units)
    h.sweepLengthInPts=synchArr(1,2)/h.nADCNumChannels;
    % the starting ticks of episodes in sample points (t0=1=beginning of
    % recording)
    h.sweepStartInPts=synchArr(:,1)*(h.synchArrTimeBase/h.fADCSampleInterval/h.nADCNumChannels);
    % recording start and stop times in seconds from midnight
    h.recTime=h.lFileStartTime+h.nFileStartMillisecs*.001;
    h.recTime=h.recTime+[0  (1e-6*(h.sweepStartInPts(end)+h.sweepLengthInPts))*h.fADCSampleInterval*h.nADCNumChannels];

  case 3
    h.dataPtsPerChan=h.lActualAcqLength/h.nADCNumChannels;
    h.dataPts=h.dataPtsPerChan*h.nADCNumChannels;
    if rem(h.dataPts,h.nADCNumChannels)>0, error('number of data points not OK'); end;
    tmp=1e-6*h.lActualAcqLength*h.fADCSampleInterval;
    if verbose
      disp('data were acquired in gap-free mode');
      disp(['total length of recording: ' num2str(tmp,'%5.1f') ' s ~ ' num2str(tmp/60,'%3.0f') ' min']);
      % 8 bytes per data point expressed in MB
      disp(['memory requirement for complete upload in matlab (gap-free data): '...
        num2str(round(8*h.lActualAcqLength/2^20)) ' MB']);
    end
    % recording start and stop times in seconds from midnight
    h.recTime=h.lFileStartTime+h.nFileStartMillisecs*.001;
    h.recTime=[h.recTime h.recTime+tmp];
  otherwise
    disp('recording mode of data must be event-driven variable-length (1), event-driven fixed-length (2) or gap-free (3) -- returning empty matrix');
    d=[];
    si=[];
end;
fclose(fid);
  
if verbose, disp(' '); end
