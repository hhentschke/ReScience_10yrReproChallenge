% this file contains
% - one partial instance of DS 
% - one partial instance of AP
% - code checking the correctness of the numerous channel indices

global DS AP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DS.abfFn=                      '02424000';    
DS.rawCh={...
'IN 0' ,1,  0,' ';...
'IN 1' ,0,  100,'pyr';...   
'IN 2' ,1,  200,' ';...     % (3)
'IN 3' ,1,  300,' ';...   
'IN 4' ,0, 400,' ';...      % (5)
'IN 7' ,0, 700,'lm';...      % (6)
'IN 12',0, 1200,'lm?';...    % (7)  
'IN 13',0, 1300,'mol?';...  % (8)    
'IN 14',0, 1400,'mol?';...  % (9)
'IN 15',1, -inf,' '...
};
DS.dpath=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AP.rawChAnNm=                DS.rawCh([3 5:9],1);
AP.rawChPrincNm=             {'IN 7'};
AP.resFn=                    [DS.abfFn '_proc_r1'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rmouse_APcheck;
rawCh=rmouse_chan;

% the following is a series of explicit comparisons

% 1. ALL channels to be analyzed
isequal(AP.chIdx,[3 5:9])
% 2. ALL LFP channels to be analyzed
isequal(AP.LFPIdx,[5:9])
isequal(AP.LFPInd,[2:6])
% 3. ALL EEG channels to be analyzed
isequal(AP.EEGIdx,[3])
isequal(AP.EEGInd,[1])
% 4. all LFP channels (2b analyzed or not)
isequal(AP.allLFPIdx,[2 5:9])
% 5. all EEG channels (2b analyzed or not)
isequal(AP.allEEGIdx,[1 3 4 10])
% 6. principal LFP channel
isequal(AP.pcIdx,[6])
isequal(AP.pcInd,[3])
% 7. indices into variables holding data from LFP computations only:
% position of analyzable LFP chans among ALL LFP chans
isequal(AP.LFPccInd,[2:6])
% ..and the left-outs
isequal(AP.LFPccOmitInd,[1]);
% position of princ chan among analyzed LFP chans
isequal(AP.LFPpcInd1,[2])
% position of princ chan among ALL LFP chans
isequal(AP.LFPpcInd2,[3])

% # of ALL channels 
isequal(AP.nAllCh,10)
% # of channels to be analyzed
isequal(AP.nCh,6)
% # of LFP channels to be analyzed
isequal(AP.nLFPCh,5)
% # of ALL LFP channels
isequal(AP.nAllLFPCh,6)
% # of EEG channels to be analyzed
isequal(AP.nEEGCh,1)
% # of ALL EEG channels
isequal(AP.nAllEEGCh,4)


% nAllCh
% 
% nCh
% AP.chIdx
% rawCh
% 
% AP.allLFPIdx
% nAllLFPCh
% 
% AP.LFPIdx
% AP.LFPInd
% nLFPCh
% 
% AP.allEEGIdx
% nAllEEGCh
% 
% AP.EEGIdx
% AP.EEGInd
% nEEGCh
% 
% AP.LFPccInd
% AP.LFPccPcInd


