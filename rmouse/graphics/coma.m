function cm=coma(type,varargin)
% ** function cm=coma(type,varargin)
% outputs custom color map, or current colormap if no legal type was
% chosen. Optional input argument 'n' specifies the number of colors
% (default 128). Output argument cm is the colormap.
% ** Note that coma may call colormap-generating functions written by
% diverse authors, notably pmkmp.m (see list of color maps below). All of
% the colormaps defined within coma.m are improvised, amateurish, playful
% colormaps without any scientific underpinning (as opposed to the concept
% underlying especially pmkmp).
%
% MONOCHROME MAPS     
%   'blues' - nomen est omen
%
% DUOTONE MAPS
%   'blackred' - nomen est omen
%   'blackgreen' - nomen est omen
%   'blackblue' - nomen est omen
%   'grayyellow' - nomen est omen
%   'blueorange' - nomen est omen
%   'turquoisamber' - shades of turquois to ...guess!
%
% MULTI-COLOR MAPS
%   'bluered' - shades of blue and red with white in between
%   'blueblackred' - shades of blue and red with black in between
%   'blueblackorange' - nomen est omen
%   'mib' - well, think of a campfire at a lake at night...
% 
% pmkmp COLOR MAPS
%   'IsoL','IsoAZ','IsoAZ180','LinearL','LinLhot','CubicYF','CubicL','Edge'
%   (see pmkmp.m for details and \plotfunc\pmkmp for further illustration)

% variable n had in a previous version been called 'ncols' - account for
% that
n=128;
pvpmod(varargin)
if exist('ncols','var')
  warning('name of optional input variable ''ncols'' was changed to ''n'' - catching this for now');
  n=ncols;
  clear ncols;
end

switch type
  case 'blackred'
    cm=(linspace(0,.3,n)')*[1 1];
    cm=[sqrt(linspace(0,1,n)') cm];
  case 'blackgreen'
    cm=(linspace(0,.3,n)')*[1 1];
    cm=[cm(:,1) sqrt(linspace(0,1,n)') cm(:,2) ];
  case 'blackblue'
    cm=(linspace(0,.3,n)')*[1 1];
    cm=[cm sqrt(linspace(0,1,n)')];
  case 'grayyellow'
    cm=(linspace(.1,1,n)')*[1 1];
    cm=[cm linspace(.1,.4,n)'];
    % cm(:,[3])=cm(:,[3]).^2;
  case 'blueorange'
    cm=(linspace(0,1,n)')*[1 1 1];    
    cm(:,1)=cm(:,1).^.8;
    cm(:,2)=sin(cm(:,2)*pi/2)*.6;
    cm(:,3)=(flipud(cm(:,3))).^2;
  case 'bluered'
    cm=(linspace(.1,1,ceil(n/2))')*[1 1 1];
    cm(:,1)=sqrt(cm(:,1));
    cm(:,[2 3])=cm(:,[2 3]).^2;
    cm=[fliplr(cm);flipud(cm)];
  case 'blueblackred'
    cm=((ceil(n/2):-1:1)')*([1 .2 .1]/ceil(n/2));
    cm(:,1)=sqrt(cm(:,1));
    % cm(:,[2 3])=cm(:,[2 3]).^2;
    cm=[fliplr(cm);flipud(cm)];
  case 'blueblackorange'
    cm=((ceil(n/2):-1:1)')*([1 .1 .1]/ceil(n/2));
    cm(:,[1 2])=sqrt(cm(:,[1 2]));
    % cm(:,[2 3])=cm(:,[2 3]).^2;
    cm=[fliplr(cm);flipud(cm)];
  case 'turquoisamber'
    cm=(linspace(.1,1,ceil(n/2))')*[1 1 1];
    cm(:,[1 2])=sqrt(cm(:,[1 2]));
    cm(:,[3])=cm(:,[3]).^2;
    cm=[fliplr(cm);flipud(cm)];
  case 'mib'
    fh=figure;
    ocm=get(gcf,'colormap');
    partCm=colormap('hot');
    % get rid of white parts
    partCm=partCm(1:round(size(partCm,1)*.9),:);
    % interpolate to get number of colors right
    [n1,n2]=size(partCm);
    for g=n2:-1:1
      cm(:,g)=interp1((1:n1)',partCm(:,g),linspace(1,n1,n/2)');
    end
    clear partCm
    % swap red and blue column, then flipud & concatenate
    cm=[flipud(cm(:,[3 2 1])); cm];
    % restore original colormap & close figure
    set(gcf,'colormap',ocm);
    close(fh)
  case 'blues'
    cm=[(linspace(0,1,n)').^2 (linspace(.2,1,n)')  (linspace(.7,1,n)')];
  case 'blues_old'
    cm=(linspace(.58,1,n)')*[.75 .9 1];
    cm(:,1)=min(cm(:,1).^2,.5);
    cm(:,2)=cm(:,2).^2;
    cm(:,3)=sqrt(cm(:,3));    
    % -------------- pmkmp colormaps -------------------------------------
  case {'IsoL','IsoAZ','IsoAZ180','LinearL','LinLhot','CubicYF','CubicL','Edge'}
    cm=pmkmp(n,type);
  otherwise
    fh=[];
    % if a colormap different from the ones defined above in is specified 
    % see whether it is a built-in one...
    try 
      ocm=colormap(type);
    catch
      % ..otherwise return the default one
      warning(['colormap ' type ' does not exist']);
      fh=figure;
      ocm=get(gcf,'colormap');
    end
    [n1,n2]=size(ocm);
    % we have to do use interpolation since all default color maps are
    % stored as arrays, not by names 
    for g=n2:-1:1
      cm(:,g)=interp1((1:n1)',ocm(:,g),linspace(1,n1,n)');
    end
    % if applicable restore original colormap & close figure
    if ~isempty(fh)
      set(fh,'colormap',ocm);
      close(fh)
    end
end
