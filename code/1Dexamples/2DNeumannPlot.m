function out = 2DNeumannPlot( V, u, v, mod_multi_idx )

out = zeros( size( u ) );
for n = 1:NSamplesPlot
    for l = 1:NSamplesPlot
        for k = 1:size( multi_idx, 1 )
            cx = multi_idx( k, 1 );
            cy = multi_idx( k, 2 );
            if( mod( cx, 2 ) && mod( cy, 2 ) )
                out( m, n, l ) = out( m, n, l ) + V( k, zeroeig( m ) ) * sin( cx * u( n, l )/2 ) .* sin( cy * v( n, l )/2 );
            elseif( mod( cx, 2 ) )
                out( m, n, l ) = out( m, n, l ) + V( k, zeroeig( m ) ) * sin( cx * u( n, l )/2 ) .* cos( cy * v( n, l )/2 ) * scaling_factor( cy + 1 );
            elseif( mod( cy, 2 ) )
                out( m, n, l ) = out( m, n, l ) + V( k, zeroeig( m ) ) * cos( cx * u( n, l )/2 ) .* sin( cy * v( n, l )/2 ) * scaling_factor( cx + 1 );
            else
                out( m, n, l ) = out( m, n, l ) + V( k, zeroeig( m ) ) * cos( cx * u( n, l )/2 ) .* cos( cy * v( n, l )/2 ) * scaling_factor( cx + 1 ) * scaling_factor( cy + 1 );
            end
        end
    end
end