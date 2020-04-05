function ticks=cont2discrete(t,resol,varargin)
% ** function ticks=cont2discrete(t,resol,varargin)
% This function accomplishes a seemingly simple job, the conversion of 
%  time expressed in real-time units (e.g. ms) to ticks given resolution 
%  'resol'. Various unexpected problems may be associated with this task 
%  - see comments at the end of this m-file for an explanation.
%
%                         >>> INPUT VARIABLES >>>
%
% NAME              TYPE/DEFAULT          DESCRIPTION
%
% t                 scalar or array       time expressed in real units
% resol             scalar, 1             sampling interval (real units)
% intv              scalar, 0             if nonzero, values in the last column of t will be 
%                                          regarded as the lower bounds of time intervals; 
%                                          the output variable 'ticks' will be adjusted 
%                                          accordingly. When converting the LENGTH of a time 
%                                          interval SET INTV TO ZERO!  
% 
%                         <<< OUTPUT VARIABLES <<<
%
% NAME             TYPE/DEFAULT          DESCRIPTION
% ticks            scalar or array       The time points/intervals expressed as ticks.
%                                         The continuous point in time t=0 will be mapped
%                                         to discrete time point ticks=1
% 
% --> see also discrete2cont, mkintrvls, regspace, embedtrc
  

intv=0;
pvpmod(varargin);
if numel(resol)~=1
  error('input parameter ''resol'' must be a scalar');
end

% the awkward construction with parentheses below is a tribute to the fact 
% that e.g.
%     floor(3000.1/0.1)
% yields 30000, which is incorrect, in contrast to either of 
%     floor(3000.1*10)
%     floor(3000.1*(1/0.1))
% which yield the correct result, 30001

ticks=floor(t*(1/resol))+1;
if intv
  if size(t,2)>2
    error('if input argument intv is nonzero input argument t must have no more than 2 columns');
  end
  ticks(:,end)=ticks(:,end)-1;  
end

% Explanation:
% The major difficulty of this task lies in the trivial but easily forgotten differences 
%  between continuous and discrete variables when it comes to the calculation of 
%  time INTERVALS. 
% Imagine you had a discretely sampled voltage trace, sampling rate 1 kHz, start
%  of recording at t=0. Question: what time intervals (expressed in real-time units)
%  do the data points represent? Simple, of course the first data point represents
%  the first millisecond, the second the second, and so on. However, if you are pressed
%  to be more exact, you will come up with time intervals of the sort [0 1[ ms, 
%  [1 2[ ms and so on: one of the interval boundaries is open (in the example above 
%  the right one, [ ). If both were closed, like in [0 1], [1 2], ... ms there would be 
%  a conceptual difficulty: e.g. time point t=1 ms would be represented by data point 1 
%  AND 2, which is impossible. Ergo: one of the interval boundaries MUST be open. 
% Now the task is to write a line of code that 
%  a) converts a point in real time (e.g. 1.3 ms) to the index of the point in the 
%  sampled voltage trace that best represents the voltage at that time. 
%  b) respects the things noted above about open interval boundaries
%  An implementation in matlab code would be
%       index=floor(time/samplingInterval)+1.
%  So, 1.3 ms would be mapped to the second point, 1.0 ms also to the second, 
%  0.9 ms to the first, and so on.
%  Now, how would we specify a time interval covering 'the first three milliseconds' 
%  of the recording in matlab code? In all likelihood, we would do it like
%              timeInterval=[0 3];
%  which is not correct, because time point 3.0 ms maps to data point 4, and 
%  what we think is an interval of 3 ms gives us back 4 data points, each of
%  which represents 1 ms worth of voltage trace. This is a principal problem
%  irrespective of the details of the implementation. The true cause is the 
%  fact that one is usually either unaware of the petty issue of the open boundary
%  and/or that sth. like
%              timeInterval=[0 3-samplingInterval/2];
%  to obtain the intended chunk of data is simply awkward. 
% Further potential complications: 
% - a commonly used way of computing discrete time from continuous time, namely
%   dividing continuous time by the sampling interval, may yield unexpected results.
%   Try, for example,
%     floor(3000.1/0.1)
%  and
%     floor(3000.1*10)
%  See any differences?
% - in large projects the conversion real time to discrete time is often 
%  done in different ways by different people. In 99% of all cases this
%  is uncritical, but in the remaining 1% this is a nasty issue that has to be dealt 
%  with.

