function out = Neumann2DPlot( V, NSamplesPlot, PlotTrig, multi_idx )

out = zeros( 1, size( PlotTrig, 2 ) );
for k = 1:size( multi_idx, 1 )
    out = out + V( k ) * PlotTrig( k, : );
end
out = reshape( out, [ NSamplesPlot NSamplesPlot ] );

% 
% for n = 1:NSamplesPlot
%     for l = 1:NSamplesPlot
%         for k = 1:size( multi_idx, 1 )
%             cx = multi_idx( k, 1 );
%             cy = multi_idx( k, 2 );
%             if( mod( cx, 2 ) && mod( cy, 2 ) )
%                 out( n, l ) = out( n, l ) + V( k ) * sin( cx * u( n, l )/2 ) .* sin( cy * v( n, l )/2 );
%             elseif( mod( cx, 2 ) )
%                 out( n, l ) = out( n, l ) + V( k ) * sin( cx * u( n, l )/2 ) .* cos( cy * v( n, l )/2 ) * scaling_factor( cy + 1 );
%             elseif( mod( cy, 2 ) )
%                 out( n, l ) = out( n, l ) + V( k ) * cos( cx * u( n, l )/2 ) .* sin( cy * v( n, l )/2 ) * scaling_factor( cx + 1 );
%             else
%                 out( n, l ) = out( n, l ) + V( k ) * cos( cx * u( n, l )/2 ) .* cos( cy * v( n, l )/2 ) * scaling_factor( cx + 1 ) * scaling_factor( cy + 1 );
%             end
%         end
%     end
% end