function [lfd,newsi,cfreq,ord]=lofi(d,si,cfreq,varargin)
% ** function [lfd,newsi,cfreq,ord]=lofi(d,si,cfreq,varargin)
%   Passes data d through butterworth lowpass filter using double reverse 
%   filtering routine 'filtfilt'. Double reverse filtering has the advantage 
%   of zero phase shift of the resulting signal. See notes at the end of 
%   this m-file for additional information
% 
%                         >>> INPUT VARIABLES >>>
% NAME              TYPE/DEFAULT         DESCRIPTION
% d                 array                data to be filtered (along columns!)
% si                scalar               sampling interval in us
% cfreq             scalar               cutoff (also termed corner, or -3 dB) frequency 
%                                         of filter in Hz
% rs                scalar, 30           minimal attenuation in stopband (dB) - see buttord
%                                         also, minimal steepness of attenuation (dB/octave) 
%                                         for frequencies between passband and stopband
% rp                scalar, 0.5          max allowed ripples in passband (dB) - see buttord
% pickf             scalar (int), 1      downsampling factor
% force_sigproc     scalar, 1            if zero, uses standard matlab routine 'filter'  
%                                         in case signal processing toolbox is not available
%                                         (a warning will be issued). If you do not want this 
%                                         behaviour, set to any nonzero value
% verbose           scalar, 1            if 1, numerical results of the calculations will be 
%                                         printed on screen. If any other nonzero value, 
%                                         in addition the filter's frequency response will 
%                                         be plotted
%
%                         <<< OUTPUT VARIABLES <<<
% NAME              TYPE/DEFAULT         DESCRIPTION
% lfd               1d- or 2d-array      filtered data 
% newsi             scalar               the sampling interval of domnsampled data in us
% cfreq             scalar               cutoff frequency in Hz resulting from computations
% ord               scalar               length of filter parameter arrays (order of filter)


% default values
rs=30;
rp=0.5;
pickf=1;
force_sigproc=1;
verbose=0;

pvpmod(varargin);
if verbose
  disp(['**** ' mfilename ':']);
end

if exist('efreq','var')
  warndlg([mfilename ': input par ''efreq '' changed to ''cfreq''']);
  creq=efreq;
end
  
% starting value for passband 'end' freq: cutoff freq
efreq=cfreq;
% rp and rs must be halved because double filtering as employed here amounts
% to doubling the effective values of rp and rs
rp=.5*rp;
rs=.5*rs;
wn=nan;

if pickf<0, error('check pick factor'); end;   
% half the original sampling rate in Hz
hfsample=.5*1e6/si;		
try
  % the -1.5 dB magnitude threshold 
  magThresh=10^(-1.5/20);
  finito=0; count=1; 
  % the factor by which deviation of current -3dB freq from targeted corner freq will be
  % multiplied and added to efreq below to make things converge
  stepFac=.4;
  % shift start and stop band ends towards lower freqs until cutoff freq is reached
  % possible future improvements: vary stepFac according to convergence behavior
  % (too slow/divergence etc.)
  while ~finito
    % edge frequency of passband normalized to hfsample, the Nyquist freq 
    wp=efreq/hfsample;		
    % stopband should start at double frequency; if, however, this is beyond
    % the Nyquist frequency, set ws to the Nyquist freq 
    ws=min(2*wp,.9999);
    [ord,wn] = buttord(wp,ws,rp,rs);
    % deviation of computed freq from -3dB freq in Hz (should initially be negative):
    fdev_3dB=cfreq-wn*hfsample;
    % generate filter parameters..
    [b,a]=butter(ord,wn);
    % ..to find the -1.5 dB freq using freqz. The -1.5 dB freq must be located between
    % cfreq and the -3 dB freq of the current run, so investigate this range
    frange=[-1:.01:1]*abs(fdev_3dB)+mean([cfreq wn*hfsample]);
    mag=abs(freqz(b,a,frange,2*hfsample));
    [nix,tmpi]=min(abs(mag-magThresh));
    % should be negative
    fdev=cfreq-frange(tmpi);
    % reduce or increase efreq accordingly by small amounts
    efreq=efreq+fdev*stepFac;    
    % minimal required accuracy: 1 %
    finito=abs(fdev/cfreq)<=.01;
    count=count+1;
  end
  cfreq=frange(tmpi);
  if verbose
    % ----- assess steepness of filter using freqz
    rs_est=[];
    % take the inner fifth of frequencies between passband and stopband
    tmp0=abs(diff([wp ws]))/5;
    tmp1=mean([wp ws])-tmp0;
    frange=([0:.1:1]*tmp0*2+tmp1)*hfsample;
    mag=abs(freqz(b,a,frange,2*hfsample));    
    % steepness of magnitude slope (in double logarithmic coordinate system)
    for i=1:length(frange)-1
      rs_est(i)=log10(mag(i)/mag(i+1))/log2(frange(i+1)/frange(i));
    end
    % multiply by 20 because that's the definition of dB
    rs_est=20*mean(rs_est);

    disp(['* Butterworth filter order: ' int2str(ord)]);
    disp(['* Iterations needed: ' int2str(count)]);    
    disp('* Filter characteristics (single filtering):')
    disp(['   Passband end/-1.5 dB/-3 dB frequencies: ' num2str(efreq,'%4.1f') '/' num2str(cfreq,'%4.1f') '/' num2str(wn*hfsample,'%4.1f') ' Hz']);
    disp(['   Mininmal atten. steepness: ' num2str(rs) ' dB/octave']);
    disp(['   Estimated maximal atten. steepness: ' num2str(rs_est,'%4.0f') ' dB/octave']);
    disp(['   Max ripples in passband: ' num2str(rp) ' dB']);
    disp('* Filter characteristics (DOUBLE filtering):')
    disp(['   -3 dB frequency: ' num2str(cfreq,'%4.1f') ' Hz']);
    disp(['   Mininmal atten. steepness: ' num2str(2*rs) ' dB/octave']);
    disp(['   Estimated maximal atten. steepness: ' num2str(2*rs_est,'%4.0f') ' dB/octave']);
    disp(['   Max ripples in passband: ' num2str(2*rp) ' dB']);
    if verbose~=1
      % show filter response in interesting freq range
      frange=min(4*efreq,hfsample)/2047*[0:2047];
      freqz(b,a,frange,2*hfsample);      
    end
  end
  if force_sigproc
    lfd=filtfilt(b,a,d);
  else
    lfd=filter(b,a,d);
  end
catch
  if force_sigproc
    error(lasterr);
  else
    % displaying the error message is necessary in case the cause is not the unavailability 
    % of the signal processing toolbox but rather some other mishap
    warndlg(lasterr);
    % in case signal processing toolbox is not available, work around..
    [fname,pname] = uigetfile('*.*','Select file for filter parameter ''a'' (ascii)');    
    a=load([pname fname],'-ascii');
    [fname,pname] = uigetfile('*.*','Select file for filter parameter ''b'' (ascii)');    
    b=load([pname fname],'-ascii');
    lfd=filter(b,a,d);
  end
end  

newsi=si;
if pickf>1,
  p=[1:pickf:size(d,1)];
  lfd=lfd(p,:,:);
  newsi=pickf*si;
  if verbose  
    disp(['New sampling interval: ' num2str(newsi) ' µs']);
  end  
end;


% NOTE: in characterizing a filter of a certain type one usually specifies the corner
% (or cutoff) frequency and the attenuation e.g. as dB per octave. Matlab offers
% functions which allow specification of either, but, unfortunately, not both at 
% the same time. This is accomplished here. The function 'buttord' expects wp and ws,
% those 'ends' of the passband and the stopband which face each other. (wp and ws are 
% termed 'corner' frequencies in the help, which is very confusing, because corner
% frequencies in the strict sense are those frequencies at which the amplitude of the
% signal is attenuated by 3 dB. Even if rp, the maximal loss in the passband, is set
% to 3, the resulting cutoff frequency (wn) does not correspond to wp). buttord
% also needs rs, the attenuation in the stopband. buttord guarantees that the filter 
% it produces has at least this attenuation in the stopband. One trick employed here 
% is the following: if we set ws to double the frequency of wp (e.g. wp=20 -> ws=40 Hz)
% then rs will also be the minimal steepness of the generated filter expressed in a 
% commonly used way: as dB per octave (octave= doubling of frequency). 
% What about the cutoff frequency? buttord calculates it on the basis of the input 
% parameters, notably the extent of the passband, but the reverse - generating a 
% filter with a certain cutoff frequency and adjusting the passband edge frequency - 
% it cannot do. This is done in the while loop above, by sequentially adjusting 
% the value of wp until wn corresponds to the desired cutoff frequency (with a 
% certain accuracy). 

