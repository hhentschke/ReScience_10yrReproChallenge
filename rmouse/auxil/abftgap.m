function dt=abftgap(abffn1,abffn2,varargin)
% ** function dt=abftgap(abffn1,abffn2,varargin)
% puts out time difference in seconds of the beginnings of two abf files
%                         >>> INPUT VARIABLES >>>
%
% NAME              TYPE/DEFAULT          DESCRIPTION
% abffn1, abffn2    char arrays
% f                 char array, 'std'     'std' - standard display
%                                         'full' - extended display
%                         <<< OUTPUT VARIABLES <<<
%
% NAME             TYPE/DEFAULT           DESCRIPTION
% dt               scalar                 time difference in s

f='std';
pvpmod(varargin);

[nix,nix2,a1]=abfload(abffn1,'info');
[nix,nix2,a2]=abfload(abffn2,'info');
% a1=abfinfo(abffn1);
% a2=abfinfo(abffn2);

if strcmpi(f,'full')
  disp('rec times file 1:')
  a1.recTime
  disp('rec times file 2:')
  a2.recTime
end
dt=a2.recTime(1)-a1.recTime(1);
disp(['time difference is ' num2str(dt,'%6.3f') ' s']);
