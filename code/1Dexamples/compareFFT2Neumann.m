function compareFFT2Neumann( idx, boo  )

% in the complex exponential basis the indices correspond to the power of
% the exponential, for the trigonometric basis the even ones correspond to
% sin and the odd ones correspond to the cosine one.

% parameters
NBases = 30;
NSamples = 1001;

% generate multi-idx to make clear what each element in row is
multi_idx = (-NBases:NBases)';
base_freqs = (-NBases ):1/2:( NBases );

A = zeros( length( multi_idx ), length( base_freqs ) );

for k = 1:length( base_freqs )
    for j = 1:length( multi_idx )
        cidx = base_freqs( k );
        tidx = multi_idx( j )/2;
        if( cidx == 0 && tidx == 0 )
            A( j, k ) = 2 * pi;
        elseif( mod( multi_idx( j ), 2 ) && ( cidx == tidx || cidx == -tidx ) )
            A( j, k ) = 1i * pi * sign( cidx * tidx );
        elseif( mod( multi_idx( j ), 2 ) )
            A( j, k ) = ( 2 * 1i * cidx * sin( pi * tidx ) * cos( pi * cidx ) )/( tidx^2 - cidx^2 );
        elseif( cidx == tidx || cidx == -tidx )
            A( j, k ) = pi;
        else
            A( j, k ) = ( -2 * cidx * cos( pi * tidx ) * sin( pi * cidx ) )/( tidx^2 - cidx^2 );
        end
    end
end

% tic
x = linspace( -pi, pi , NSamples );
outFF = zeros( NSamples, 1 );
outF = zeros( NSamples, 1 );
outC = zeros( NSamples, 1 );

% h = A * f1vec';
h = A * boo;
for n = 1:NSamples
    outF( n ) = ( 2 * pi )^-1 * sum( boo .* exp( 1i * x( n ) * base_freqs( : ) ) );
    outFF( n ) = ( 2 * pi )^-1 * sum( boo( mod( base_freqs, 1 ) == 0 ) .* exp( 1i * x( n ) * base_freqs( mod( base_freqs, 1 ) == 0 )' ) );
    for m = 1:length( multi_idx )
        if ( mod( multi_idx( m ), 2 ) )
            outC( n ) = outC( n ) + h( m ) * sin( multi_idx( m ) * x( n )/2 );
        else
            outC( n ) = outC( n ) + h( m ) * cos( multi_idx( m ) * x( n )/2 );
        end
    end
    outC( n ) = outC( n ) * ( 2 * pi )^-2;
end
%  x = linspace( -pi, pi , NSamples );
% outF = zeros( NSamples, 1 );
% for n = 1:NSamples
%     outF( n ) = ( 2 * pi )^-1 * sum( boo .* exp( 1i * x( n ) * sampled_freqs' ) );
% end
%
idx
if ( mod( idx, 2 ) )
    figure; hold on;
%     blah = sin( x ) .* ( idx/2 ) .* cos( idx * x/2 ) + 1/2 .* sin( idx * x/2 ) .* cos( x );
    blah = ( 1 - x ) .* ( pi^2 - x.^2 ) .* ( idx/2 ) .* cos( idx * x/2 ) + 1/2 .* sin( idx * x/2 ) .* ( x .* ( 3 .* x - 2 ) - pi^2 );
    plot( x, blah, 'LineWidth', 1 );
    plot( x, real( outF ), 'r', 'LineWidth', 1 );
    plot( x, real( outFF ), 'c', 'LineWidth', 1 );
    plot( x, real( outC ), 'g' );
else
    figure; hold on;
%     blah = sin( x ) .* ( -idx/2 ) .* sin( idx * x/2 ) + 1/2 .* cos( idx * x/2 ) .* cos( x );
    blah = ( 1 - x ) .* ( pi^2 - x.^2 ) .* ( -idx/2 ) .* sin( idx * x/2 ) + 1/2 .* cos( idx * x/2 ) .* ( x .* ( 3 .* x - 2 ) - pi^2 );
    plot( x, blah, 'LineWidth', 1 );
    plot( x, real( outF ), 'r', 'LineWidth', 1 );
    plot( x, real( outFF ), 'c', 'LineWidth', 1 );
    plot( x, real( outC ), 'g' );
end


max( blah' - real( outF ) )
max( real( outF ) - real( outC ) )

