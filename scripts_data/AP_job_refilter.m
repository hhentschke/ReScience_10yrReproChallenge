% this is a partial AP (Analysis Parameter) set to be used in batch processing 


AP.job={...
  'gen_btsl',...            % generate behavioral time stamp list
  'filter&hilbert',...      % generate the various streams (theta, gamma, etc.)
};

% save plot(s) to file in what format? '-dpsc2' (color postscript level 2) is
% strongly recommended. multiple formats can be specified; in this case they must 
% be listed in a cell array, like {'-dpsc2','-djpeg90'}. set to [] if plots 
% shall not be saved.
AP.printas=                  [];{'-djpeg90'};

