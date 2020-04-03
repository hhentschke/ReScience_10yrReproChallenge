% this is a partial AP (Analysis Parameter) set to be used in batch processing 

AP.job={...
  'det_artifacts',...       % detect artifacts in neural signals & generate time stamp list
  'gen_btsl',...            % generate behavioral time stamp list
  'sumFig',...              % produce summary plots
  'tcFig',...               % plots time course of segment-wise computed, selected parameters 
  'tcPrincComp',...         % computes & plots principal components from tupels of segment-wise computed, selected parameters   
};

AP.job={...
  'det_artifacts',...       % detect artifacts in neural signals & generate time stamp list
  'gen_btsl',...            % generate behavioral time stamp list
  'sumFig',...              % produce summary plots
  'tcFig',...               % plots time course of segment-wise computed, selected parameters 
};

AP.job={...
  'det_artifacts',...       % detect artifacts in neural signals & generate time stamp list
  'gen_btsl',...            % generate behavioral time stamp list
  'sumFig',...              % produce summary plots
};

% save plot(s) to file in what format? '-dpsc2' (color postscript level 2) is
% strongly recommended. multiple formats can be specified; in this case they must 
% be listed in a cell array, like {'-dpsc2','-djpeg90'}. set to [] if plots 
% shall not be saved.
AP.printas=                  {'-djpeg90'};[];

