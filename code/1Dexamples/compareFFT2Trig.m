% in the complex exponential basis the indices correspond to the power of
% the exponential, for the trigonometric basis the even ones correspond to
% sin and the odd ones correspond to teh cosine one.

A = zeros( ( 2 * NBases + 1 ) );

for k = 1:length( multi_idx )
    for j = 1:length( multi_idx )
        cidx = multi_idx( k );
        tidx = multi_idx( j )/2;
        if ( mod( multi_idx( j ), 2 ) )
            A( j, k ) = ( 2 * tidx * sin( ( pi * tidx ) ) * cos( pi * cidx ) )/( tidx^2 - cidx^2 );
        elseif( cidx == tidx )
            A( j, k ) = 1i * pi * sign( cidx * tidx );
        elseif( cidx == -tidx )
            A( j, k ) = 1i * pi * sign( cidx * tidx );
        end
    end
end

% tic
x = linspace( -pi, pi , NSamples );
outF = zeros( NSamples, 1 );
outC = zeros( NSamples, 1 );
h = A * f1vec';
for n = 1:NSamples
    outF( n ) = ( 2 * pi )^-1 * sum( f1vec' .* exp( 1i * x( n ) * multi_idx( : ) ) );
    for m = 1:length( multi_idx )
        if ( mod( multi_idx( m ), 2 ) )
            outC( n ) = outC( n ) + h( m ) * cos( multi_idx( m ) * x( n )/2 );
        else
            outC( n ) = outC( n ) + h( m ) * sin( multi_idx( m ) * x( n )/2 );
        end
    end
    outC( n ) = outC( n ) * ( 2 * pi )^-1;
end

