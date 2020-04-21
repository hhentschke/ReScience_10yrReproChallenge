function [p,F,radj1,radj2]=curvecomp(crf1,crf2,crf_c,df1,df2,varargin)
% ** function [p,F,radj1,radj2]=curvecomp(crf1,crf2,crf_c,df1,df2,varargin)
% tests difference between two dose-response curves for statistical
% significance. Fits of a model to the data are required. Reference: H.
% Motulsky, Analyzing data with GraphPad Prism (Version 3), p. 231 ff (This
% is a software manual freely available in the Internet). The same
% information, but scattered about chapters, can be found in H. Motulsky &
% A. Christopoulos, Fitting Models to biological data using linear and
% nonlinear regression, chaps 22 and 27.
% hh Aug 03, modified June 04 
%                         >>> INPUT VARIABLES >>>
%
% NAME         TYPE/DEFAULT       DESCRIPTION
% crf1, crf2   column arrays      crf = 'concentration', 'response', 'fit'
%                                 the two concentration response curves to be compared
%                                 they must hold the following numbers:
%                                 column 1: concentrations
%                                 column 2: response (e.g., average firing rate at these concentrations)
%                                 column 3: fit (e.g. values of fitted Hill curve at these concentrations)
% crf_c        column array       same as crf1 and crf2 but obtained from the COMBINATION of both 
%                                 (procedure: combine the experimental data from both data sets and 
%                                 fit the same curve type to this merged data set)
% df1, df2     scalars            degrees of freedom for each of the fitted curves
%                                 calculation: number of data points - number of fitted variables
%                                 in the case of e.g. a Hill fit the number of fitted variables 
%                                 is 3
%                             
%
%                         <<< OUTPUT VARIABLES <<<
%
% NAME          TYPE/DEFAULT           DESCRIPTION
% F             scalar                 F-value
% p             scalar                 p-value
% radj1, radj2  scalars                'degrees of freedom-adjusted R square', a measure of 
%                                      goodness of fit of each of the individual data sets 
%                                      (the closer it is to 1, the better)
% 

% ----- default values & varargin -----
verbose=0;
pvpmod(varargin);
if verbose, disp(['**** ' mfilename ':']); end;


% simple checks
[rows1,cols1]=size(crf1);
[rows2,cols2]=size(crf2);
[rows_c,cols_c]=size(crf_c);
if cols1~=3 | cols2 ~= 3 | cols_c ~= 3
  error('at least one of the input variables crf1, crf2, crf_c has more or less than 3 columns');
end
if rows1<2 | rows2<2 | rows_c<2
  error('at least one of the input variables crf1, crf2, crf_c has less than 2 rows (=data points)');
end
if rows_c ~= rows1+rows2
  error('crf_c MUST have as many rows as crf1 and crf2 togehter');
end

if df1>=rows1 | df2>=rows2
  error('one or both of the data sets (crf) has less data points (=rows) than degrees of freedom (df), which is impossible');
end
% the number of variables in the fit = number of data points - number of degrees of freedom
nVar1=rows1-df1;
nVar2=rows2-df2;
if nVar1 ~= nVar2
  error('number of degrees of freedom in crf1 and crf2 does not match');
end


% 0. for both data sets, compute a handful of terms a) as descriptors of the goodness of fit 
% and b) for comparison with combined model further below
% + SSE, summed square of residuals
sse1=sum((crf1(:,2)-crf1(:,3)).^2);
sse2=sum((crf2(:,2)-crf2(:,3)).^2);
% + Root Mean Squared Error (also called standard error of the regression)
rmse1=sqrt(sse1/df1);
rmse2=sqrt(sse2/df2);
% + total sum of squares
sst1=sum((crf1(:,2)-mean(crf1(:,2))).^2);
sst2=sum((crf2(:,2)-mean(crf2(:,2))).^2);
% + R square: 
r1=1-sse1/sst1;
r2=1-sse2/sst2;
% + degrees of freedom-adjusted R square, also termed 'adjusted coefficient 
% of determination', in German 'korrigiertes Bestimmtheitsmass':
% ** note: values below differ from the output produced by the matlab fit
% function (statistics toolbox V 5.1): within the fit function, radj is
% computed as  1-sse*(row-1)/(sst*df)
radj1=1-sse1*(rows1-1)/(sst1*(df1-1));
radj2=1-sse2*(rows2-1)/(sst2*(df2-1));


% 1. add individual SSE
ss_summed=sse1+sse2;
% ..and degrees of freedom
df_summed=df1+df2;

% 2. compute same for combined data set
ss_combo=sum((crf_c(:,2)-crf_c(:,3)).^2);
df_combo=df1+df2+nVar1;

% 3. compute F- and p-values
F=((ss_combo-ss_summed)/(df_combo-df_summed))/(ss_summed/df_summed);
p=1-fcdf(F,df_combo-df_summed,df_summed);
