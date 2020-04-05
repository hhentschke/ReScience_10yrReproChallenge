function d=i16load(file,varargin)
% ** function d=i16load(file,varargin)
% Reads file, interpreting contents as binary 2 byte signed integers
%
%                         >>> INPUT VARIABLES >>>
%
% NAME             TYPE/DEFAULT          DESCRIPTION
% file             string                file name
% nPts             scalar, inf           number of points to read starting with 
%                                         first point (inf=all)
%             ---- OR -----
% intv             2-element array       first and last points to read
% byteOrd          string, 'ieee-le'     byte order
%                                        'ieee-le' - little-endian (win)
%                                        'ieee-be' - big-endian (unix, labview)
% verbose          scalar,1              if nonzero, details of loading will be
%                                         printed on screen
%
%                         <<< OUTPUT VARIABLES <<<
%
% NAME             TYPE/DEFAULT          DESCRIPTION
% d                column array          the data

% default values
verbose=1;
% byte order
byteOrd='ieee-le';
nPts=inf;
pvpmod(varargin);
if verbose, disp(['** ' mfilename ':']); end;

% 2 byte/data point
bytePP=2;
% the last position the file pointer can move to according to user's specification
if exist('intv','var')
  nPts=diff(intv)+1;
else
  intv=[1 nPts];
end

if verbose, disp(['opening file ' file '..']); end;
fid=fopen(file,'r',byteOrd);
if intv(2)~=inf
  if fseek(fid,bytePP*intv(2),'bof')~=0,		
    message=ferror(fid);
    error(['requested points do not match file size (' int2str(nPts) ');(' message ')']);
  end;
end
% rewind to intv(1)
fseek(fid,(intv(1)-1)*bytePP,'bof');
d=fread(fid,nPts,'int16');
fclose(fid);