function rmouse_ini
% This file contains some default parameters (individual preferences) and 
% 'working prameters' (WP) that are specific to the machine the analysis is 
% run on

global WP

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%  GENERAL
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
set(0,'units','pixels')
% The progress of computations is by default printed on screen and in 
% a logfile. If WP.verbose is set to 1, much more detail will be given.
% Don't set to any value other than 0 or 1.
WP.verbose=0;
% which version of matlab?
WP.mver=ver;
WP.mver=WP.mver(strmatch('matlab',lower({WP.mver.Name}),'exact')).Version;
% ..and is java loaded?
WP.javaEnabled=version('-java');
if ~isempty(strfind(WP.javaEnabled,'not enabled'))
  WP.javaEnabled=0;
else
  WP.javaEnabled=1;
end
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%  PERFORMANCE & MEMORY ISSUES
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% The maximal amount of memory (in MegaBytes) the entirety of large variables 
% may occupy at any time during analysis. If set too low, some of the analysis 
% jobs must switch to a different (slower) working mode (involving many read 
% and write operations to/from the original and intermediate results files). 
% On the other hand, if this value is too high, the OS will spend much time
% swapping memory to disk or run out of memory. A rule of thumb-value is two
% thirds of the computer's physical RAM. For further orientation, the 
% amount of memory needed to hold a specific piece of raw data in the matlab 
% workspace is
% [number of channels] * [60 s/min * data length (min) * sampling rate (Hz)]
%  * [8 (bytes per data point (the standard matlab double))]
% So, assuming a sampling rate of 1000 Hz, 
% number of channels      length of data   memory requirement
% 1                       30 min           ~ 14 MB
% 8                       30 min           ~110 MB
% 16                      30 min           ~220 MB
% 
% Now assume that for a recording of 30 min with 16 channels the major results 
% variable r will take up about 150 Mb of RAM. So, on a machine with 512 Mb 
% of RAM it would be advisable to set WP.maxRAM to a little more than 
% 220 + 150 Mb because that will allow matlab to hold raw data of all channels
% in memory at once, which speeds things up.

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%  PATHS
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% root directories for data may differ from machine to machine. If DS.dpath does not
% contain a drive (c: or d:), then WP.rootPath will be pre-pended to DS.dpath

if exist('c:\grndhh','file'), 
  WP.maxRAM=1200;
  WP.rootPath='d:\rmouse';  
elseif exist('c:\smlhh','file'), 
  WP.maxRAM=80;
  WP.rootPath='e:\rmouse';  
elseif exist('c:\lv_dual','file'), 
  WP.maxRAM=750;
  WP.rootPath='x:\rmouse';
elseif strcmpi(getenv('computername'),'hh64')
  WP.maxRAM=2000;
  WP.rootPath='d:\rmouse';
elseif strcmpi(getenv('computername'),'HH-I7')
  WP.maxRAM=20000;
  WP.rootPath='d:\_data\rmouse\ten-years';
else 
  error('root path is undefined (see rmouse_ini.m)');  
end

