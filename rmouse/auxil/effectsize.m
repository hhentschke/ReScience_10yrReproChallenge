function [es,varargout]=effectsize(x,y,varargin)
% ** function [es,varargout]=effectsize(x,y,varargin)
% computes two measures of effect size, Cohen's d or Hedges' d, between two
% samples.
% All input parameters except x and y are optional and must be specified as
% parameter/value pairs, e.g. as in
%          effectsize(x,y,'type','cohend');
%
%                         >>> INPUT VARIABLES >>>
%
% NAME             TYPE/DEFAULT          DESCRIPTION
% x, y           arrays                data samples (observations in rows;
%                                         in case of several columns
%                                         comparisons will be made among
%                                         matching columns)
% type             char, 'hedgesd'       'hedgesd' or 'cohend'
% pair             scalar, 0             if data are paired (repeated
%                                         measures) set to nonzero
% nBoot            scalar, nan           if set to positive value: number
%                                         of bootstrap repetitions for
%                                         computation of CI
%                         <<< OUTPUT VARIABLES <<<
%
% NAME             TYPE/DEFAULT           DESCRIPTION
% es                                      effect size
% ci               varargout{1}           95% confidence intervals
%
% NOTE: this code is outdated. Better use function mes.m.



% ----- default values & varargin -----
type='hedgesd';
pair=0;
nBoot=nan;
confLevel=.05;
pvpmod(varargin);

% ----- a few 'constants':
% minimal number of bootstrapping runs below which a warning will be issued
% (and bottstrapping will be skipped) - set to a value of your taste
minNBootstrap=100;

% ----- a few input checks
if ~ismember(type,{'cohend','hedgesd'})
  error('illegal type of effect size specified');
end

if strcmp(type,'cohend') & pair
  error('sorry, Cohen''s d not implemented for paired data - please use Hedges'' d');
end

[nRowX nColX]=size(x);
[nRowY nColY]=size(y);
if nColX~=nColY
  error('input variables x and y must have the same number of columns');
end
if pair
  if nRowX~=nRowY
    error('for paired data input variables x and y must have the same number of rows');
  end
end

% check output
if nargout>1
  doComputeCi=true;
else
  doComputeCi=false;
end

% check bootstrapping settings
doBoot=false;
if doComputeCi
  if isfinite(nBoot)
    if nBoot>=minNBootstrap;
      doBoot=true;
    else
      warning('number of bootstrap repetitions is not adequate - not bootstrapping');
    end
  end
end


if doBoot
  % *** bootstrap
  % retain original data in separate variables
  x_orig=x;
  y_orig=y;
  % preallocate es and ci in case of multiple comparisons
  es=repmat(nan,1,nColX);
  ci=repmat(nan,2,nColX);
  h=waitbar(0,'Bootstrapping...');
  % loop over number of comparisons
  for g=1:nColX
    % pick corresponding columns, omit nans
    x=x_orig(isfinite(x_orig(:,g)),g);
    y=y_orig(isfinite(y_orig(:,g)),g);
    % generate bootstrap indices
    [bootstat,bootsam]=bootstrp(nBoot,[],x);
    % first column of x: original data
    x=[x x(bootsam)];
    % §§§ if data are paired use same bootsample for y; if they are unpaired
    % compute new one
    if ~pair
      [bootstat,bootsam]=bootstrp(nBoot,[],y);
    end
    % first column of y: original data
    y=[y y(bootsam)];
    % preparatory works: compute n, means, variances
    [n1,m1,s1]=prepEsComp(x);
    [n2,m2,s2]=prepEsComp(y);
    
    % compute effect size from original and bootstrapped data
    switch type
      case 'cohend'
        cur_es=cohend(n1,n2,m1,m2,s1,s2);
      case 'hedgesd'
        cur_es=hedgesg(x,y,n1,n2,m1,m2,s1,s2,pair);
    end
    % compute confidence intervals
    ci(1:2,g)=bootconfi(cur_es,confLevel);
    % retain first element of cur_es; this is the actual es
    es(g)=cur_es(1);
    waitbar(g/nColX,h);
  end
  % kill waitbar window
  delete(h);
  varargout{1}=ci;
  
else
  % *** no bootstrapping
  % preparatory works: compute n, means, variances
  [n1,m1,s1]=prepEsComp(x);
  [n2,m2,s2]=prepEsComp(y);
  % for the computation of confidence intervals below use inverse
  % cumulative t distribution with the following degrees of freedom:
  if pair
    % paired: n=n1-1=n2-1
    ciFac= -tinv(confLevel/2,n1-1);
  else
    % independent: n=[n1+n2-2]
    ciFac= -tinv(confLevel/2,n1+n2-2);
  end
  
  switch type
    case 'cohend'
      es=cohend(n1,n2,m1,m2,s1,s2);
      if doComputeCi
        varargout{1}=cohend_ci(n1,n2,es,ciFac);
      end
    case 'hedgesd'
      [es,r]=hedgesg(x,y,n1,n2,m1,m2,s1,s2,pair);
      if doComputeCi
        varargout{1}=hedgesg_ci(n1,n2,es,ciFac,r,pair);
      end
  end
end


% ========================= LOCAL FUNCTIONS ===============================
function [n,m,s]=prepEsComp(x)
% preparatory computations: n, means, variances, etc.
n=sum(isfinite(x),1);
m=nanmean(x,1);
s=nanvar(x,0,1); % var and std are by default divided by n-1

% Cohen's d
% *** Note: all functions below are able to deal with ARRAYS (as opposed to
% scalars) as inputs. Thus, they are able to deal, in a vectorized manner,
% with multiple comparisons (i.e. if input arguments x and y into
% effectsize.m are multi-column arrays) as well as with bootstrapping.
function es=cohend(n1,n2,m1,m2,s1,s2)
es=(m1-m2)./sqrt(((n1-1).*s1 + (n2-1).*s2)./(n1+n2));

% confidence intervals for Cohen's d, analytical (Nakagawa & Cuthill 2007)
function ci=cohend_ci(n1,n2,es,ciFac)
se=sqrt( (n1+n2-1)./(n1+n2-3).*(4./(n1+n2)).*(1+es.^2/8));
ci=cat(1,es-ciFac.*se,es+ciFac.*se);

% Hedges' g plus t statistics which are needed for analytical confidence
% intervals (below)
function [es,r]=hedgesg(x,y,n1,n2,m1,m2,s1,s2,pair)
if pair
  % paired data
  [h,p,tci,tst]=ttest(x,y);
  % compute pairwise correlation coeffs 
  r=repmat(nan,1,size(x,2));
  for g=1:size(x,2)
    r(g)=corr(x(:,g),y(:,g));
  end
  % note that n=n1=n2
  es=tst.tstat.*sqrt((2-2*r)./n1);
else
  r=nan;
  es=(m1-m2)./sqrt(((n1-1).*s1 + (n2-1).*s2)./(n1+n2-2));
end
% correct for bias
es=es.*(1-(3./(4*n1+4*n2-9)));

% confidence intervals for Hedges' g, analytical (Nakagawa & Cuthill 2007)
function ci=hedgesg_ci(n1,n2,es,ciFac,r,pair)
if pair
  % paired data - note that n=n1=n2
  se=sqrt((2-2*r)./n1 + es.^2./(2*n1-1));
else
  se=sqrt((n1+n2)./(n1.*n2) + (es.^2./(2*n1+2*n2-4)));
end
ci=cat(1,es-ciFac.*se,es+ciFac.*se);

% compute confidence intervals from arrays of effect size measures
% generated from bootstrapped data
function ci=bootconfi(es,confLevel)
ci=prctile(es,[confLevel/2  1-confLevel/2]'*100);


