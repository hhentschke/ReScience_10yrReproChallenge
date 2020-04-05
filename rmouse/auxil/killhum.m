function d=killhum(d,si,f)
% **function d=killhum(d,si,f)
% Passes data d through butterworth bandstop filter with approximate corner
% frequencies [f-.5 f+.5] Hz using double reverse filtering routine
% 'filtfilt'. si is the sampling interval of the data in us.

warning('THIS FUNCTION IS A DINOSAUR - consider using more adequate functions (e.g. elim_hum)')

rp=.5;			% max attenuation in passband, don't change
rs=90;			% attenuation in dB per octave

hfsample=(1e6/si)/2;		% this must be HALF the original sampling rate IN HZ
fcut=[f-0.5 f+0.5];
wp=fcut/hfsample;			
ws(1,1)=0.5*wp(1,1);
ws(1,2)=2*wp(1,2);		
[n,wn] = buttord(wp,ws,rp,rs);
[b,a]=butter(n,wn,'stop');
% plot frequency response?
% freqz(b,a,512,2*hfsample);
disp(['Butterworth natural 3dB frequencies: ' num2str(wn*hfsample) ' Hz']);
disp(['Order ' num2str(n)]);
disp(['Attenuation: ' num2str(rs) ' dB/octave']);
disp(['Ripples: ' num2str(rp) ' dB in passband']);

d=filtfilt(b,a,d);