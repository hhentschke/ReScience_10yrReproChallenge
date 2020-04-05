function etslconst
% This file contains constants needed when dealing with variables of 
% 'extended time stamp list' (etsl) format. An etsl is a N by 4 array
% containing events as defined by a time stamp, duration, amplitude and type.
% --> see primer_evltsl.rtf for additional explanation

% ** changed global to assignin nov 2010 **

% global etslc;


% the columns:
% - time stamps
etslc.tsCol=1;
% - duration
etslc.durCol=2;
% - amplitude
etslc.amplCol=3;
% - tags
etslc.tagCol=4;
% the number of columns expected
etslc.nCol=4;

% tags attached to events in tagCol;
% events which are not O.K. in terms of data integrity 
% (particularly missing raw data preventing cutout generation etc.)
etslc.corrupt=-1;
% events which are not OK in terms of data quality (artifacts etc.)
etslc.reject=-2;
% not (yet) applicable
etslc.na=NaN; 

assignin('caller','etslc',etslc);