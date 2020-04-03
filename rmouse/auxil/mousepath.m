function mousePath
% ** function mousePath
% opens a simple questdlg, asks for the principal data path to be used for
% analysis of 'rmouse' data, sets it to the chosen one and removes the
% other

% set path to which data directory?
ddir={'beta3_wtko','WTb3N265M'};
curWorkDir=questdlg('Choose principal path','mouse question',ddir{1},ddir{2},ddir{2});

if exist('c:\grndhh','file'), 
  mouseRootP='c:\rmouse';  
elseif exist('c:\smlhh','file'), 
  mouseRootP='d:\rmouse';  
elseif exist('c:\dual_chore','file'), 
  mouseRootP='d:\rmouse';  
elseif exist('c:\lv_dual','file'), 
  mouseRootP='x:\rmouse';  
elseif exist('c:\matlabor','file'), 
  mouseRootP='d:\rmouse';  
elseif exist('c:\dagobert','file'), 
  mouseRootP='g:\rmouse';  
else 
  errordlg('unknown machine, no path set or removed');
end;

addpath([mouseRootP '\' curWorkDir]);
path2remove=[mouseRootP '\'  ddir{setdiff([1 2],strmatch(curWorkDir,ddir))}];
if ~isempty(strfind(path,[path2remove pathsep]))
  rmpath(path2remove)
end
cd([mouseRootP '\' curWorkDir]);
