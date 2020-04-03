function rmouse_APcheck

global DS AP

% here, Analysis Structure AP is checked for integrity (missing and superfluous fields)
tmpAP=AP;
AP=[];
anpar_template;
s=fieldnames(AP);
s2=fieldnames(tmpAP);

% 1. set all fields of template AP to empty matrix
for fi=1:length(s)
  eval(['c=iscell(AP.' s{fi} ');']);
  if c, eval(['AP.' s{fi} '={[]};']);
  else eval(['AP.' s{fi} '=[];']);
  end
end
templateAP=AP;

% 2. identify and create missing fields
[c,ix]=setdiff(s,s2);
if ~isempty(c)
  disp('appending fields to AP:');
  disp(c);
  for fi=1:length(ix)
    eval(['tmpAP.' c{fi} '=AP.' c{fi} ';']);
  end
end
% 3. eliminate fields not needed
s=fieldnames(templateAP);
s2=fieldnames(tmpAP);
[c,ix]=setdiff(s2,s);
if ~isempty(c)
  warning('deleting outdated fields from AP:');
  disp(c);
  tmpAP=rmfield(tmpAP,c);
end
% finally, order fields
AP=orderfields(tmpAP,templateAP);