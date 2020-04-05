function t=discrete2cont(ticks,resol,varargin)
% ** function t=discrete2cont(ticks,resol,varargin)
% This function accomplishes a simple job, the conversion of time expressed 
%  in ticks (also termed bins, or 'sampled points') to real-time units (e.g. ms) given 
%  resolution 'resol'. 
% This function is the counterpart to cont2discrete, see comments there.
% 
%                         >>> INPUT VARIABLES >>>
%
% NAME              TYPE/DEFAULT          DESCRIPTION
%
% ticks             scalar or array       time expressed in ticks
% resol             scalar, 1             sampling interval (real units)
% intv              scalar, 0             if nonzero, values in the first column will be 
%                                         regarded as the upper bounds of time intervals; 
%                                         the output variable 'ticks' will be adjusted 
%                                         accordingly
%
%                         <<< OUTPUT VARIABLES <<<
%
% NAME             TYPE/DEFAULT          DESCRIPTION
% ticks             scalar or array       the time points/intervals expressed as ticks.
%                                          The discrete data point with index 1 is mapped
%                                          to the continuous point in time t=0
%
% --> see also cont2discrete, embedtrc, mkintrvls, regspace


intv=0;
pvpmod(varargin);

if numel(resol)~=1
  error('input parameter ''resol'' must be a scalar');
end
if ~isequaln(round(ticks),ticks)
  error('ticks must be integer numbers');
end

t=(ticks-1)*resol;
if intv
  if size(ticks,2)>2
    error('if input argument intv is nonzero input argument ticks must have no more than 2 columns');
  end
  t(:,end)=t(:,end)+resol;
end
