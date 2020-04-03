% produces all megaplots one could ever wish for
opt=3;

% pn - qualifier for file names
switch opt
  case 0
    pn='foo';
    wavefm=0;
    % figure directory
    fidi='\beta3_wtko\figures';
  case 1
    pn='WTKO_line';
    wavefm=1;
    % figure directory
    fidi='\beta3_wtko\figures';
  case 2
    pn='WT_iState';
    pn='KO_iState';
    wavefm=0;
    % figure directory
    fidi='\beta3_wtko\figures';
  case 3
    pn='WT';
    pn='KO';    
    pn='WT_Atrop_Ctrl';        

    pn='WT_Atrop_App';                    
    
    masterAP='AP_beta3_wtko';
    wavefm=0;
    % figure directory
    fidi='\beta3_wtko\figures';
  case 4
    pn='WT';    
    pn='KI';
    wavefm=0;    
    fidi='\WTb3N265M\figures';
  case 5
    pn='WT_ctrl_15';    
    pn='WT_eto_15';    

    pn='KI_ctrl_15';
    pn='KI_eto_15';

    wavefm=0;    
    fidi='\WTb3N265M\figures';
    masterAP='AP_beta3_wtki';
end


behav={'immobile','exploring'};


% all gamma-related
rv={'gaCCMn','gaeCCMn','thgaeCCMn',...
    'gaCCPeakMn','gaeCCPeakMn',...
    'gaCCPeakTMn','gaeCCPeakTMn'};
% all
rv={'deCCMn','thCCMn','gaCCMn','thLoeCCMn','thHieCCMn','gaeCCMn','detheCCMn','thgaeCCMn',...
    'deCCPeakMn','thCCPeakMn','gaCCPeakMn','thLoeCCPeakMn','thHieCCPeakMn','gaeCCPeakMn',...
    'deCCPeakTMn','thCCPeakTMn','gaCCPeakTMn','thLoeCCPeakTMn','thHieCCPeakTMn','gaeCCPeakTMn',...
    'rawCohMn','fComod','fComodP',...
    'rawDePEMn','rawThPEMn','rawThNarrowPEMn','rawBePEMn','rawGaPEMn','rawRiPEMn'};

rv={'thCCMn','gaCCMn','gaeCCMn','thgaeCCMn',...
    'thCCPeakMn','gaCCPeakMn','gaeCCPeakMn',...
    'thCCPeakTMn','gaCCPeakTMn','gaeCCPeakTMn',...
    'rawThPEMn','rawThNarrowPEMn','rawBePEMn','rawGaPEMn'};

% all new
rv={'rawgaeCohPeak','rawgaeCohPeakF','rawgaeCohTh'};

r_megaplot(masterAP,'behav',behav,'rv',rv,'wavefm',wavefm,'pn',pn,'figdir',fidi);
