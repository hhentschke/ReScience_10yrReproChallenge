% collect ANPARS etc.
% substance should be either F6 or Iso or F6iso
substance='Iso';
%substance='F6iso';
% conc must be any combination of 'control' , 'drug' , 'wash'
 conc={'Control','Drug'};
conc={'Control','Drug','Wash'};

% ------------------------------------------------
ANPAR=[];DSET=[];AP=[];DS=[];A=[];D=[];
eval(['collect_' substance ';']);
A=ANPAR;
D=DSET;
mAP='AP_ratF6iso';
% figure directory
fidi='c:\data\Rat invivo Data\F6invivo\figures';


behav={'immobile','exploring'};
rv={'deCCMn','thCCMn','gaCCMn','thLoeCCMn','thHieCCMn','gaeCCMn','detheCCMn','thgaeCCMn',...
    'deCCPeakMn','thCCPeakMn','gaCCPeakMn','thLoeCCPeakMn','thHieCCPeakMn','gaeCCPeakMn',...
    'deCCPeakTMn','thCCPeakTMn','gaCCPeakTMn','thLoeCCPeakTMn','thHieCCPeakTMn','gaeCCPeakTMn'};

for bi=1:length(behav)
  for rvi=1:length(rv)
    for ci=1:length(conc)
      switch lower(conc{ci})
        case 'control'
          % reduce to control measurements
          DSET=D(1,:);
          ANPAR=A(1,:);
        case 'drug'
          % reduce to F6 application
          DSET=D(2,:);
          ANPAR=A(2,:);
        case 'wash'
          % reduce to wash measurements
          DSET=D(3,:);
          ANPAR=A(3,:);
        otherwise
          error('check spelling')
      end
      r_megaplot(mAP,'behav',behav(bi),'rv',rv(rvi),'pn',[substance conc{ci}],'figdir',fidi);
    end
  end
end
