function c=strmwrite(fid,d)
% ** function c=strmwrite(fid,d)
% writes data (streams) as scaled 16 bit integers 
global DS
c=fwrite(fid,round((d/DS.nsRng)*2^15),'int16');
