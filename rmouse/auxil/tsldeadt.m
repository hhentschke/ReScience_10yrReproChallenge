function [ntsl,ix]=tsldeadt(tsl,deadT,varargin)
% ** function [ntsl,ix]=tsldeadt(tsl,deadT,varargin)
%   crops a time stamp list such that each time stamp has a minimum 
%   distance in time (dead time) to the previous time stamp. Within each 
%   group of crowded time stamps, the first will survive, the others will
%   go.
%
%        ** time unit must be the same for all variables **
%
%                    >>> INPUT VARIABLES >>>
%
% NAME          TYPE/DEFAULT         DESCRIPTION
%
% tsl           column array         time stamp list
% deadT         scalar               the dead time
% include_1st   scalar, 1            if nonzero, first event in list 
%                                    will be included               
%
%                    <<< OUTPUT VARIABLES <<<
%
% NAME          TYPE/DEFAULT          DESCRIPTION
%
% ntsl          column array          the new, cropped time stamp list
% ix            column array          logical index into tsl to elements 
%                                     that made it into ntsl

% hh 08/03, modified 1/07

include_1st=1;
pvpmod(varargin)

if isempty(tsl),
  error(['tsl is empty']);
else
  if deadT<0
    error('negative dead time');
  elseif deadT==0
    disp([mfilename ': dead time is zero - no cropping']);
    ntsl=tsl;
    ix=(1:numel(tsl))';
  else
    ix=[logical(include_1st); diff(tsl)>=deadT];
    ntsl=tsl(ix);
  end
end
