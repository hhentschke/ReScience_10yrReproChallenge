function [funH, fitParNames, beta, ds1ix, ds2ix, ds12ix]=fit2rmousepar(ds1, ds2, ds12, rv)
% fit2rmousepar has been added in April 2020 for the then years
% reproducibility challenge, the major difference being that it works
% without the curve fitting toolbox.

switch rv

  
  case 'rawGaCentroidMn'
    % including data for x=0 makes sense for this CC parameter
    ds1ix=find(ds1(:,1)>-inf);
    ds2ix=find(ds2(:,1)>-inf);
    ds12ix=find(ds12(:,1)>-inf);
    % third-order polynomial
    ft_ = fittype('a + b*x + c*(x*x) + d*(x*x*x)' ,...
      'dependent',{'y'},'independent',{'x'},...
      'coefficients',{'a','b','c','d'});
    fo_ = fitoptions('method','NonlinearLeastSquares');
    % starting values for parameters
    st_ = [60 0 0 0];

      
  case {'thNegPeakCvAMn','thNegPeakCvIPIMn','thPosPeakCvAMn','thPosPeakCvIPIMn','gaePosPeakCvIPIMn'}
    % including data for x=0 makes sense for this parameter
    ds1ix=find(ds1(:,1)>-inf);
    ds2ix=find(ds2(:,1)>-inf);
    ds12ix=find(ds12(:,1)>-inf);
    if 1
      % gaussian + offset (parameter b1 could be fixed at ~40 without any
      % discernible loss in the quality of the fit!)
      ft_ = fittype('a1*exp(-b1*(x-c1)^2) + o' ,...
        'dependent',{'y'},'independent',{'x'},...
        'coefficients',{'a1','b1','c1','o'});
      fo_ = fitoptions('method','NonlinearLeastSquares');
      % starting values for parameters
      st_ = [.2 30 .4 .2];
    else
      % second-order polynomial
      ft_ = fittype('a + b*x + c*(x*x)' ,...
        'dependent',{'y'},'independent',{'x'},...
        'coefficients',{'a','b','c'});
      fo_ = fitoptions('method','NonlinearLeastSquares');
      % starting values for parameters
      st_ = [.5 .5 -1];
    end
    
  case 'gaePosPeakCvAMn'
    % including data for x=0 makes sense for this CC parameter
    ds1ix=find(ds1(:,1)>-inf);
    ds2ix=find(ds2(:,1)>-inf);
    ds12ix=find(ds12(:,1)>-inf);
    % second-order polynomial
    ft_ = fittype('a + b*x + c*(x*x)' ,...
      'dependent',{'y'},'independent',{'x'},...
      'coefficients',{'a','b','c'});
    fo_ = fitoptions('method','NonlinearLeastSquares');
    % starting values for parameters
    st_ = [.35 -.2 .3];



  case {'thCCPosPeakDecayMn'}
    % including data for x=0 makes sense for this CC parameter
    ds1ix=find(ds1(:,1)>-inf);
    ds2ix=find(ds2(:,1)>-inf);
    ds12ix=find(ds12(:,1)>-inf);
    fo_ = fitoptions('method','NonlinearLeastSquares');
    if 1
      % second-order polynomial
      ft_ = fittype('a + b*x + c*(x*x)' ,...
        'dependent',{'y'},'independent',{'x'},...
        'coefficients',{'a','b','c'});
      % starting values for parameters
      st_ = [.3 -15 15];
    else
      % third-order polynomial
      ft_ = fittype('a + b*x + c*(x*x) + d*(x*x*x)' ,...
        'dependent',{'y'},'independent',{'x'},...
        'coefficients',{'a','b','c','d'});
      % starting values for parameters
      st_ = [.3 -1 1 1];
    end
    
  
  case {'rawDePEMn','rawPMnDeP'}
    % including data for x=0 makes sense for power measurements
    ds1ix=find(ds1(:,1)>-inf);
    ds2ix=find(ds2(:,1)>-inf);
    ds12ix=find(ds12(:,1)>-inf);
    if 1
      % exponential + offset
      % ** the Radj for this fit is so bad that a linear function seems more
      % appropriate
      ft_ = fittype('a*exp(-b*x) + o' ,...
        'dependent',{'y'},'independent',{'x'},...
        'coefficients',{'a','b','o'});
      % starting values for parameters - potentially sensitive!
      st_ = [.1 10 .005];
    else
      ft_ = fittype('a*x+b' ,...
        'dependent',{'y'},'independent',{'x'},...
        'coefficients',{'a','b'});
      st_ = [0 .1];
    end
    fo_ = fitoptions('method','NonlinearLeastSquares');

  case {'rawThPEMn','rawThNarrowPEMn','rawPMnThP','rawPMnThNarrowP'}
    % including data for x=0 makes sense for power measurements
    ds1ix=find(ds1(:,1)>-inf);
    ds2ix=find(ds2(:,1)>-inf);
    ds12ix=find(ds12(:,1)>-inf);
    % new: exponential + offset
    ft_ = fittype('a*exp(-b*x) + o' ,...
      'dependent',{'y'},'independent',{'x'},...
      'coefficients',{'a','b','o'});
    fo_ = fitoptions('method','NonlinearLeastSquares');
    % starting values for parameters - potentially sensitive!
    st_ = [.2 10 .005];

  case {'rawBePEMn','rawBeNarrowPEMn'}
    % including data for x=0 makes sense for power measurements
    ds1ix=find(ds1(:,1)>-inf);
    ds2ix=find(ds2(:,1)>-inf);
    ds12ix=find(ds12(:,1)>-inf);
    % new: exponential + offset
    ft_ = fittype('a*exp(-b*x) + o' ,...
      'dependent',{'y'},'independent',{'x'},...
      'coefficients',{'a','b','o'});
    fo_ = fitoptions('method','NonlinearLeastSquares');
    % starting values for parameters - potentially sensitive!
    st_ = [.05 10 .005];

  case {'rawGaPEMn','rawGaNarrowPEMn','gaeThPEMn'}
    % including data for x=0 makes sense for power measurements
    ds1ix=find(ds1(:,1)>-inf);
    ds2ix=find(ds2(:,1)>-inf);
    ds12ix=find(ds12(:,1)>-inf);
    if 0
      % old: Gaussian type + offset -> quite good if only dorsal sites are fitted
      ft_ = fittype('a*exp(-b*(x)^2) + o' ,...
        'dependent',{'y'},'independent',{'x'},...
        'coefficients',{'a','b','o'});
      fo_ = fitoptions('method','NonlinearLeastSquares');
      % starting values for parameters - potentially sensitive!
      st_ = [.02 10 .005];
    else
      % new: exponential + offset
      ft_ = fittype('a*exp(-b*x) + o' ,...
        'dependent',{'y'},'independent',{'x'},...
        'coefficients',{'a','b','o'});
      fo_ = fitoptions('method','NonlinearLeastSquares');
      % starting values for parameters - potentially sensitive!
      st_ = [.02 10 .005];
    end
  
  case {'gaeCCPeakMn'}
    % including data for x=0 does not make sense for CC parameters
    ds1ix=find(ds1(:,1)~=0);
    ds2ix=find(ds2(:,1)~=0);
    ds12ix=find(ds12(:,1)~=0);
      % gaussian * (slowing cos + offset), force first hump to reside at x=0 with amplitude of 1
      ft_ = fittype('exp(-b1*(x)^2) * (a2*cos(b2*x^c2)+(1-a2))' ,...
        'dependent',{'y'},'independent',{'x'},...
        'coefficients',{'b1','a2','b2','c2'});
      fo_ = fitoptions('method','NonlinearLeastSquares');
      % starting values for parameters - it is of utmost importance that b2, the cosine
      % freq, has a good starting value
      st_ = [1 .4 12 .9];
      
  case {'thgaeCCPeakTStd','thgaeCCPeakPhaseStd'}
    % including data for x=0 makes sense for thgaeCC
    ds1ix=find(ds1(:,1)>-inf);
    ds2ix=find(ds2(:,1)>-inf);
    ds12ix=find(ds12(:,1)>-inf);
    if 1
      % logistic + parabola - works quite well!
      ft_ = fittype('a/(1+exp(-b*(x+c)))+d*x^2' ,...
        'dependent',{'y'},'independent',{'x'},...
        'coefficients',{'a','b','c','d'});
      % starting values for parameters
      st_ = [20 20 -.3 -20];
    else
      % gaussian + linear
      ft_ = fittype('a1*exp(-b1*(x-c1)^2) + a2*x' ,...
        'dependent',{'y'},'independent',{'x'},...
        'coefficients',{'a1','b1','c1','a2'});
      % starting values for parameters
      st_ = [10 20 .3 0];
    end
    fo_ = fitoptions('method','NonlinearLeastSquares');

  case {'thCCPeakTStd','thCCPeakPhaseStd','gaeCCPeakTStd','gaeCCPeakPhaseStd'}
    % including data for x=0 does not make sense for CC parameters
    ds1ix=find(ds1(:,1)~=0);
    ds2ix=find(ds2(:,1)~=0);
    ds12ix=find(ds12(:,1)~=0);
    if 1
      % logistic + parabola - works quite well!
      ft_ = fittype('a/(1+exp(-b*(x+c)))+d*x^2' ,...
        'dependent',{'y'},'independent',{'x'},...
        'coefficients',{'a','b','c','d'});
      % starting values for parameters
      st_ = [20 20 -.3 -20];
    else
      % gaussian + linear
      ft_ = fittype('a1*exp(-b1*(x-c1)^2) + a2*x' ,...
        'dependent',{'y'},'independent',{'x'},...
        'coefficients',{'a1','b1','c1','a2'});
      st_ = [10 20 .3 5];
    end
    fo_ = fitoptions('method','NonlinearLeastSquares');
   
    
  case {'gaCCPeakMn'}
    % including data for x=0 does not make sense for CC parameters
    ds1ix=find(ds1(:,1)~=0);
    ds2ix=find(ds2(:,1)~=0);
    ds12ix=find(ds12(:,1)~=0);
    % Gaussian + 3rd order polynomial, force hump to reside at x=0 with amplitude of 1
    ft_ = fittype('a1*exp(-b1*(x)^2) + a2*x^2 + a3*x^3 + 1-a1' ,...
      'dependent',{'y'},'independent',{'x'},...
      'coefficients',{'a1','b1','a2','a3'});
    fo_ = fitoptions('method','NonlinearLeastSquares');
    % fo_ = fitoptions('method','NonlinearLeastSquares',...
    % 'Lower',[.0001 .0001 -50 -50],'Upper',[1000 100 50 100]);
    % starting values for parameters
    st_ = [.9 10 .2 -1];
    
  case {'thCCPeakMn'}
    % including data for x=0 does not make sense for CC parameters
    ds1ix=find(ds1(:,1)~=0);
    ds2ix=find(ds2(:,1)~=0);
    ds12ix=find(ds12(:,1)~=0);
    if 0
      % Gaussian + 3rd order polynomial, force hump to reside at x=0 with amplitude of 1
      ft_ = fittype('a1*exp(-b1*(x)^2) + a2*x^2 + a3*x^3 + 1-a1' ,...
        'dependent',{'y'},'independent',{'x'},...
        'coefficients',{'a1','b1','a2','a3'});
      % starting values for parameters
      st_ = [2 10 80 -60];
      % 'Lower',[.0001 .0001 -50 -50],'Upper',[1000 100 50 100]);
      fo_ = fitoptions('method','NonlinearLeastSquares');
    elseif 0
    % third-order polynomial, force y=1 at x=0
    ft_ = fittype('1 + b*x + c*(x*x) + d*(x*x*x)' ,...
      'dependent',{'y'},'independent',{'x'},...
      'coefficients',{'b','c','d'});
      fo_ = fitoptions('method','NonlinearLeastSquares');
      % starting values for parameters
      st_ = [0 0 0];
    else
      % Gaussian + 2nd order polynomial, force y=1 at x=0
      ft_ = fittype('a1*exp(-b1*(x)^2) + a2*x^2 + 1-a1' ,...
        'dependent',{'y'},'independent',{'x'},...
        'coefficients',{'a1','b1','a2'});
      % starting values for parameters
      st_ = [1 10 2 ];
      fo_ = fitoptions('method','NonlinearLeastSquares');

    end

  case {'thgaeCCPeakMn','thgaeCCEnvPeakMn','rawgaeCohPeak'}
    % including data for x=0 makes sense for this CC parameter
    ds1ix=find(ds1(:,1)>-inf);
    ds2ix=find(ds2(:,1)>-inf);
    ds12ix=find(ds12(:,1)>-inf);
    
    funH = @(par, x) par(1) + par(2)*x + par(3)*x.^2;
    % parameter names
    fitParNames = ["poly0", "poly1", "poly2"];
    % starting values for parameters
    beta = [0.15, -.3, 1];
    
  case {'thgaeCCZScore'}
    % including data for x=0 makes sense for this CC parameter
    ds1ix=find(ds1(:,1)>-inf);
    ds2ix=find(ds2(:,1)>-inf);
    ds12ix=find(ds12(:,1)>-inf);
    if 0
      % OLD: sum of two gaussians + offset
      ft_ = fittype('a1*exp(-b1*(x)^2) + a2*exp(-b2*(x-c2)^2)' ,...
        'dependent',{'y'},'independent',{'x'},...
        'coefficients',{'a1','b1','a2','b2','c2'});
      fo_ = fitoptions('method','NonlinearLeastSquares');
      % starting values for parameters
      st_ = [.2 40 .15 4 .6 ];
    else
      % second-order polynomial
      ft_ = fittype('a + b*x + c*(x*x)' ,...
        'dependent',{'y'},'independent',{'x'},...
        'coefficients',{'a','b','c'});
      fo_ = fitoptions('method','NonlinearLeastSquares');
      % starting values for parameters
      st_ = [4 -15 15];
    end

  case {'gaeCCPeakTMn','gaeCCPeakPhaseMn'}
    % including data for x=0 does not make sense for CC parameters
    ds1ix=find(ds1(:,1)~=0);
    ds2ix=find(ds2(:,1)~=0);
    ds12ix=find(ds12(:,1)~=0);
    %
    ft_ = fittype('a1*exp(-b1*(x-.3)^2) + a2*x' ,...
      'dependent',{'y'},'independent',{'x'},...
      'coefficients',{'a1','b1','a2'});
    fo_ = fitoptions('method','NonlinearLeastSquares');
    % fo_ = fitoptions('method','NonlinearLeastSquares',...
    % 'Lower',[.0001 .0001 -50 -50],'Upper',[1000 100 50 100]);
    % starting values for parameters
    if strcmpi(rv,'gaeCCPeakTMn')
      st_ = [5 50 1];
    else
      st_ = [5/15 50 1/15];
    end
      
    
  case {'gaCCPeakTMn'}
    % including data for x=0 does not make sense for CC parameters
    ds1ix=find(ds1(:,1)~=0);
    ds2ix=find(ds2(:,1)~=0);
    ds12ix=find(ds12(:,1)~=0);
    % logistic type  (problem: f(0) ~=0, f'(0) ~=0)
    ft_ = fittype('a/(1 + c*exp(-b*(x)))' ,...
      'dependent',{'y'},'independent',{'x'},...
      'coefficients',{'a','b','c'});
    fo_ = fitoptions('method','NonlinearLeastSquares');,...
      % starting values for parameters - potentially sensitive!
    st_ = [10 10 30];

  case {'thCCPeakTMn','thCCPeakPhaseMn'}
    % including data for x=0 does not make sense for CC parameters
    ds1ix=find(ds1(:,1)~=0);
    ds2ix=find(ds2(:,1)~=0);
    ds12ix=find(ds12(:,1)~=0);
    % logistic type  (problem: f(0) ~=0, f'(0) ~=0)
    ft_ = fittype('a/(1 + c*exp(-b*(x)))' ,...
      'dependent',{'y'},'independent',{'x'},...
      'coefficients',{'a','b','c'});
    fo_ = fitoptions('method','NonlinearLeastSquares');,...
      % starting values for parameters - potentially sensitive!
    if strcmpi(rv,'thCCPeakTMn')
      st_ = [40 10 30];
    else
      st_ = [40/15 10 30];
    end

  case {'thgaeCCPeakTMn','thgaeCCPeakPhaseMn'}
    % including data for x=0 does make sense for this particular CC!
    ds1ix=find(ds1(:,1)>=-inf);
    ds2ix=find(ds2(:,1)>=-inf);
    ds12ix=find(ds12(:,1)>=-inf);
    if 1
      % logistic type with offset and shifted on x-axis + xxx
      % to account for the 'dip' ca. 300 um dorsal of slm
      ft_ = fittype('a/(1 + exp(-b*(x+d)))+o +c*x' ,...
        'dependent',{'y'},'independent',{'x'},...
        'coefficients',{'a','b','d','c','o'});
      % fo_ = fitoptions('method','NonlinearLeastSquares',...
      %  'Lower',[.0001 .0001 -50 -50],'Upper',[1000 100 50 100]);
      fo_ = fitoptions('method','NonlinearLeastSquares');
      % starting values for parameters - potentially sensitive!
      if strcmpi(rv,'thgaeCCPeakTMn')
        st_ = [50 25 -.4  -5 -10];
      else
        st_ = [50/15 25 -.4  -5/15 -10/15];
      end
    else
      % logistic type with offset and shifted on x-axis + xxx
      % to account for the 'dip' ca. 300 um dorsal of slm
      ft_ = fittype('a/(1+exp(-b*(x+d)))+c*x' ,...
        'dependent',{'y'},'independent',{'x'},...
        'coefficients',{'a','b','d','c'});
      % fo_ = fitoptions('method','NonlinearLeastSquares',...
      %  'Lower',[.0001 .0001 -50 -50],'Upper',[1000 100 50 100]);
      fo_ = fitoptions('method','NonlinearLeastSquares');
      % starting values for parameters - potentially sensitive!
      if strcmpi(rv,'thgaeCCPeakTMn')
        st_ = [50 25 -.4  -5];
      else
        st_ = [50/15 25 -.4  -5/15 ];
      end
    end

  case 'rawPPeakMn'
    % including data for x=0 makes sense for power measurements
    ds1ix=find(ds1(:,1)>-inf);
    ds2ix=find(ds2(:,1)>-inf);
    ds12ix=find(ds12(:,1)>-inf);
    % Gaussian type + offset -> quite good if only dorsal sites are fitted
    ft_ = fittype('a*exp(-b*(x)^2) + o' ,...
      'dependent',{'y'},'independent',{'x'},...
      'coefficients',{'a','b','o'});
    fo_ = fitoptions('method','NonlinearLeastSquares');
    % starting values for parameters - potentially sensitive!
    st_ = [.02 10 .005];

  case 'rawPPeakTMn'
    % including data for x=0 makes sense for power measurements
    ds1ix=find(ds1(:,1)>-inf);
    ds2ix=find(ds2(:,1)>-inf);
    ds12ix=find(ds12(:,1)>-inf);
    %
    ft_ = fittype('a*x+b' ,...
      'dependent',{'y'},'independent',{'x'},...
      'coefficients',{'a','b'});
    fo_ = fitoptions('method','NonlinearLeastSquares');
    % starting values for parameters - potentially sensitive!
    st_ = [0 8];

  case 'rawPMnPeak'
    % including data for x=0 makes sense for power measurements
    ds1ix=find(ds1(:,1)>-inf);
    ds2ix=find(ds2(:,1)>-inf);
    ds12ix=find(ds12(:,1)>-inf);
    % Gaussian type + offset -> quite good if only dorsal sites are fitted
    ft_ = fittype('a*exp(-b*(x)^2) + o' ,...
      'dependent',{'y'},'independent',{'x'},...
      'coefficients',{'a','b','o'});
    fo_ = fitoptions('method','NonlinearLeastSquares');
    % starting values for parameters - potentially sensitive!
    st_ = [.02 10 .005];

  case 'rawPMnPeakT'
    % including data for x=0 makes sense for power measurements
    ds1ix=find(ds1(:,1)>-inf);
    ds2ix=find(ds2(:,1)>-inf);
    ds12ix=find(ds12(:,1)>-inf);
    %
    ft_ = fittype('a*x+b' ,...
      'dependent',{'y'},'independent',{'x'},...
      'coefficients',{'a','b'});
    fo_ = fitoptions('method','NonlinearLeastSquares');
    % starting values for parameters - potentially sensitive!
    st_ = [0 8];

  case {'rawCohMnDe'}  % *** not yet scrutinzed
    % including data for x=0 does not make sense for CC parameters
    ds1ix=find(ds1(:,1)~=0);
    ds2ix=find(ds2(:,1)~=0);
    ds12ix=find(ds12(:,1)~=0);
    % Gaussian + 3rd order polynomial, force hump to reside at x=0 with amplitude of 1
    ft_ = fittype('a1*exp(-b1*(x)^2) + a2*x^2 + a3*x^3 + 1-a1' ,...
      'dependent',{'y'},'independent',{'x'},...
      'coefficients',{'a1','b1','a2','a3'});
    fo_ = fitoptions('method','NonlinearLeastSquares');
    % fo_ = fitoptions('method','NonlinearLeastSquares',...
    % 'Lower',[.0001 .0001 -50 -50],'Upper',[1000 100 50 100]);
    % starting values for parameters
    st_ = [2 10 80 -60];

  case {'rawCohMnTh','rawCohMnThNarrow','gaeCohMnThNarrow'}
    % including data for x=0 does not make sense for CC parameters
    ds1ix=find(ds1(:,1)~=0);
    ds2ix=find(ds2(:,1)~=0);
    ds12ix=find(ds12(:,1)~=0);
    if 0
      % OLD, but not bad: gaussian * (cos + offset), force first hump to reside at x=0 with amplitude of 1
      ft_ = fittype('exp(-b1*(x)^2) * (a2*cos(b2*x^c2)+(1-a2))' ,...
        'dependent',{'y'},'independent',{'x'},...
        'coefficients',{'b1','a2','b2','c2'});
      fo_ = fitoptions('method','NonlinearLeastSquares');
      % fo_ = fitoptions('method','NonlinearLeastSquares',...
      % 'Lower',[.0001 .0001 -50 -50],'Upper',[1000 100 50 100]);
      % starting values for parameters - it is of utmost importance that b2, the cosine
      % freq, has a good starting value
      st_ = [1 .2 7.8 1];
    else
      % Gaussian + 3rd order polynomial, force hump to reside at x=0 with amplitude of 1
      ft_ = fittype('a1*exp(-b1*(x)^2) + a2*x^2 + a3*x^3 + 1-a1' ,...
        'dependent',{'y'},'independent',{'x'},...
        'coefficients',{'a1','b1','a2','a3'});
      % starting values for parameters
      st_ = [2 10 80 -60];
      % 'Lower',[.0001 .0001 -50 -50],'Upper',[1000 100 50 100]);
      fo_ = fitoptions('method','NonlinearLeastSquares');
    end

  case {'rawCohMnGa'}
    % including data for x=0 does not make sense for CC parameters
    ds1ix=find(ds1(:,1)~=0);
    ds2ix=find(ds2(:,1)~=0);
    ds12ix=find(ds12(:,1)~=0);
    % Gaussian + 3rd order polynomial, force hump to reside at x=0 with amplitude of 1
    ft_ = fittype('a1*exp(-b1*(x)^2) + a2*x^2 + a3*x^3 + 1-a1' ,...
      'dependent',{'y'},'independent',{'x'},...
      'coefficients',{'a1','b1','a2','a3'});
    fo_ = fitoptions('method','NonlinearLeastSquares');
    % fo_ = fitoptions('method','NonlinearLeastSquares',...
    % 'Lower',[.0001 .0001 -50 -50],'Upper',[1000 100 50 100]);
    % starting values for parameters
    st_ = [2 10 80 -60];

  case {'rawCohMnRi'}  % *** not yet scrutinzed
    % including data for x=0 does not make sense for CC parameters
    ds1ix=find(ds1(:,1)~=0);
    ds2ix=find(ds2(:,1)~=0);
    ds12ix=find(ds12(:,1)~=0);
    % Gaussian + 3rd order polynomial, force hump to reside at x=0 with amplitude of 1
    ft_ = fittype('a1*exp(-b1*(x)^2) + a2*x^2 + a3*x^3 + 1-a1' ,...
      'dependent',{'y'},'independent',{'x'},...
      'coefficients',{'a1','b1','a2','a3'});
    fo_ = fitoptions('method','NonlinearLeastSquares');
    % fo_ = fitoptions('method','NonlinearLeastSquares',...
    % 'Lower',[.0001 .0001 -50 -50],'Upper',[1000 100 50 100]);
    % starting values for parameters
    st_ = [2 10 80 -60];

  
  otherwise
    % all other cases: linear so as not to stop program
    % include data for x=0
    ds1ix=find(ds1(:,1)>-inf);
    ds2ix=find(ds2(:,1)>-inf);
    ds12ix=find(ds12(:,1)>-inf);
    %
    ft_ = fittype('a*x+b' ,...
      'dependent',{'y'},'independent',{'x'},...
      'coefficients',{'a','b'});
    fo_ = fitoptions('method','NonlinearLeastSquares');
    % starting values for parameters - potentially sensitive!
    st_ = [0 2];
    warning(['no fit specially designed for ' rv ' - employing linear fit']);
end % switch

