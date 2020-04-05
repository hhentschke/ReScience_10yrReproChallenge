function fh=mkfig(name)

fh=[];
n=[];
sz = 'b';
% pvpmod(varargin)  

if isempty(fh) || ~isgraphics(fh)
  if isempty(n) || ~isfinite(n)
    fh=figure;
  else
    fh=figure(n);
  end
else
  figure(fh);
end

scs=get(0,'screensize');
marg=round(scs(4)/40);
switch sz
  case 'b'
    set(fh,'position',[scs(1)+marg  floor(scs(4)*.25)-marg   scs(3)*.6              floor(scs(4)*.75)-2*marg]);
  case 'h'
    set(fh,'position',[scs(1)+marg  floor(scs(4)*.45)-marg   scs(3)*1-2*marg        floor(scs(4)*.5)]);
  case 'v'
    set(fh,'position',[scs(1)+marg  scs(2)+2*marg            floor(scs(3)*.5)-marg  scs(4)-5*marg]);
  case 'max'
    set(fh,'position',[scs(1)+marg  scs(2)+2*marg            scs(3)*1-2*marg        scs(4)-5*marg]);
  case 'min'
    set(fh,'position',[scs(1)+marg  floor(scs(4)*.8)-marg    floor(scs(3)*.2)-marg  floor(scs(4)*.2)-marg]);
end

if ~isempty(name)
  set(fh,'name',name);
end