function rawCh=rmouse_chan

global DS AP WP

% This file generates about every kind of supplementary variable/index to
% channels one would ever want to have

% There are two basic types of indices to channels, '...Idx' and '...Ind'. 
% The difference: 
%   - '...Idx' always refer to DS.rawCh, that is, the set of channels
%   listed in DS (regardless of their type)
%   - '...Ind' refer to AP.rawChAnNm, that is, to the set of channels to be
%   analyzed. AP.rawChAnNm itself depends on (=is a subset of) DS.rawCh,
%   which means that all ...Ind variables and derivates MUST NOT be used as
%   indices to DS.rawCh 

% There is one more type:
%   - '...ccInd' refer to AP.LFPIdx (all LFP channels to be analyzed).
%   rmouse is designed predominantly for the analysis of LFP data;
%   therefore, index arrays are needed which pick out analyzable LFP
%   channels among all LFP channels 

% # of ALL channels 
AP.nAllCh=size(DS.rawCh,1);
% # of channels to be analyzed (this excludes the scoring and wheel velocity channels)
AP.nCh=length(AP.rawChAnNm);

% GLOBAL (=related to DS.rawCh) and/or LOCAL indices to 
% 1. ALL channels to be analyzed
AP.chIdx=[];
% 2. ALL LFP channels to be analyzed
AP.LFPIdx=[];
AP.LFPInd=[];
% 3. ALL EEG channels to be analyzed
AP.EEGIdx=[];
AP.EEGInd=[];
% 4. all LFP channels (2b analyzed or not)
AP.allLFPIdx=find([DS.rawCh{:,2}]==0);
% 5. all EEG channels (2b analyzed or not)
AP.allEEGIdx=find([DS.rawCh{:,2}]==1);
% 6. principal LFP channel
AP.pcIdx=strmatch(AP.rawChPrincNm{1},DS.rawCh(:,1),'exact');
AP.pcInd=strmatch(AP.rawChPrincNm{1},AP.rawChAnNm,'exact');
if isempty(AP.pcIdx), error('principal channel must be included in set of channels to be analyzed'); end
% % warn if principal channel is not labeled 'lm'
% if ~strcmpi(DS.rawCh{AP.pcIdx,4},'lm')
%   warnstr=[DS.abfFn ': principal channel in AP not labeled ''lm'' in DS'];
%   warning(warnstr);
%   warndlg(warnstr);
% end    
% check channel type - attention, in the current version scoring channels are not allowed!
for i=1:AP.nCh      
  AP.chIdx(i)=strmatch(AP.rawChAnNm{i},DS.rawCh(:,1),'exact');
  if DS.rawCh{AP.chIdx(i),2}==0
    AP.LFPIdx=[AP.LFPIdx AP.chIdx(i)];
    AP.LFPInd=[AP.LFPInd i];
  elseif DS.rawCh{AP.chIdx(i),2}==1
    AP.EEGIdx=[AP.EEGIdx AP.chIdx(i)];
    AP.EEGInd=[AP.EEGInd i];    
  else
    error('illegal data type specified in DS.rawCh (2nd column)')
  end
end
% 7. indices into variables holding data from LFP computations only:
% position of analyzable LFP chans among ALL LFP chans
[nix,AP.LFPccInd]=intersect(AP.allLFPIdx,AP.LFPIdx);
% ..and the inverse: the left-outs
[nix,AP.LFPccOmitInd]=setdiff(AP.allLFPIdx,AP.LFPIdx);
% position of princ chan among analyzed LFP chans
AP.LFPpcInd1=find(AP.LFPIdx==AP.pcIdx);
% position of princ chan among ALL LFP chans
AP.LFPpcInd2=find(AP.allLFPIdx==AP.pcIdx);

% # of LFP channels to be analyzed
AP.nLFPCh=length(AP.LFPIdx);
% # of ALL LFP channels
AP.nAllLFPCh=length(AP.allLFPIdx);
% # of EEG channels to be analyzed
AP.nEEGCh=length(AP.EEGIdx);
% # of ALL EEG channels
AP.nAllEEGCh=length(AP.allEEGIdx);

% linear index in 2D CC LFP cell arrays accessing main diagonal
AP.dixie=1:AP.nAllLFPCh+1:AP.nAllLFPCh^2;
% linear index in 2D CC LFP cell arrays accessing upper triangular part including main diagonal
co=combin(AP.nAllLFPCh,'autoC',1);
AP.trixie=sub2ind([AP.nAllLFPCh,AP.nAllLFPCh],co(:,1),co(:,2));
% linear index in 2D CC LFP cell arrays accessing upper triangular part EXCLUDING main diagonal
co=combin(AP.nAllLFPCh,'autoC',0);
AP.strixie=sub2ind([AP.nAllLFPCh,AP.nAllLFPCh],co(:,1),co(:,2));

% position of channel to start thgaeCC with among ALL LFP chans
if ~isempty(AP.OPT_thgaeCCStartChNm{1})
  tmpi=strmatch(AP.OPT_thgaeCCStartChNm{1},DS.rawCh(:,1),'exact');
  AP.OPT_thgaeCCStartChInd=find(AP.allLFPIdx==tmpi);
  if isempty(AP.OPT_thgaeCCStartChInd)
    error('AP.OPT_thgaeCCStartChNm is specified but does not match any of the channels to be analyzed');
  end
else
  AP.OPT_thgaeCCStartChInd=[];
end
if ~isempty(AP.OPT_thgaeCCAttractLag) && length(AP.OPT_thgaeCCAttractLag) ~= AP.nLFPCh
  error('AP.OPT_thgaeCCAttractLag has too many or too few entries');
end
% position of channel to start detheCC with among ALL LFP chans
if ~isempty(AP.OPT_detheCCStartChNm{1})
  tmpi=strmatch(AP.OPT_detheCCStartChNm{1},DS.rawCh(:,1),'exact');
  AP.OPT_detheCCStartChInd=find(AP.allLFPIdx==tmpi);
  if isempty(AP.OPT_detheCCStartChInd)
    error('AP.OPT_detheCCStartChNm is specified but does not match any of the channels to be analyzed');
  end
else
  AP.OPT_detheCCStartChInd=[];
end
if ~isempty(AP.OPT_detheCCAttractLag) && length(AP.OPT_detheCCAttractLag) ~= AP.nLFPCh
  error('AP.OPT_detheCCAttractLag has too many or too few entries');
end

% spacing of all LFP electrodes: 
% WP.elx is the true depth in mm (mm look better on plots than um) 
% WP.felx is the 'functional' depth, i.e. the true depth linearly scaled
% (such that distance PYR-LM =1); elx and felx refer to the same channels
% WP.xt* are the labels for the various x axes
WP.elx=(cat(1,DS.rawCh{AP.allLFPIdx,3})-DS.rawCh{AP.pcIdx,3})*.001;
% get rid of annoying rounding errors
WP.elx=0.01*round(WP.elx*100);
WP.xtl1=repmat(' ',[AP.nAllLFPCh 4]);
WP.xtl2=repmat(' ',[AP.nAllLFPCh 2]);
xtlspace=floor(AP.nAllLFPCh/8)+1;
for i=1:xtlspace:AP.nAllLFPCh, 
  WP.xtl1(i,:)=sprintf('%+1.1f',WP.elx(i)); 
  % extract numerical part of channel name
  tmpChnm=DS.rawCh{AP.allLFPIdx(i),1};
  WP.xtl2(i,:)=sprintf('%2i',str2double(tmpChnm(~isletter(tmpChnm))));         
end;
WP.felx=nan*ones(AP.nAllLFPCh,1);
ix1=strmatch('pyr',DS.rawCh(:,4),'exact');
ix2=strmatch('lm',DS.rawCh(:,4),'exact');
if length(ix1)==1 && length(ix2)==1
  WP.felx=(cat(1,DS.rawCh{AP.allLFPIdx,3})-DS.rawCh{ix2,3})/(DS.rawCh{ix2,3}-DS.rawCh{ix1,3});  
elseif length(ix1)~=1 && length(ix2)==1
  % skip this annoying warning until functional depth axis is really needed
  % warning('''pyr'' not specified or redundant in DS.rawCh - cannot set up functional depth axis');
  % ** presence of 'lm' has already been checked above
end  

% ---- define all streams that rmouse can deal with
WP.strm={...
  'sansDelta',...
  'delta',...
  'theta',...
  'thetaEnv',...
  'thetaLo',...
  'thetaLoEnv',...
  'thetaHi',...
  'thetaHiEnv',...
  'beta',...
  'gamma',...
  'gammaEnv',...
  'gammaNarrow',...
  'gammaNarrowEnv',...
  'ripple',...
  'rippleEnv'};

% ---- make sure all of the streams listed in AP.strm are legal 
if ~isempty(setdiff(AP.strm,WP.strm))
  error('illegal stream in AP.strm (see rmouse_chan.m for a list of legal streams)');
end

% generate struct rawCh which holds - for all channels 2b analyzed -
% channel name, file names of preprocessed streams of these channels and
% channel type
for i=1:AP.nCh
  rawCh(i).nm=AP.rawChAnNm{i};
  switch DS.rawCh{AP.chIdx(i),2}
    case 0
      rawCh(i).dType='LFP';  
    case 1
      rawCh(i).dType='EEG';  
  end
  tmpNm=[DS.abfFn '_' AP.rawChAnNm{i}];
  for g=1:numel(WP.strm)
    rawCh=setfield(rawCh,{i},[WP.strm{g} 'Fn'],[tmpNm '_' WP.strm{g} '.i16']);
  end
end

% name of file for theta peaks and gammaEnv peaks
WP.thetaPFn=[AP.strmDir '\' DS.abfFn '_thPeaks.mat'];
WP.gammaEnvPFn=[AP.strmDir '\' DS.abfFn '_gaePeaks.mat'];

