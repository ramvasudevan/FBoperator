function out = Neumann2DPlot( V, NSamplesPlot, PlotTrig )

out = V' * PlotTrig;
out = reshape( out, [ NSamplesPlot NSamplesPlot ] );


% for k = 1:size( multi_idx, 1 )
%     out = out + V( k ) * PlotTrig( k, : );
% end