function compareFFT2Trig( idx, boo  )

% in the complex exponential basis the indices correspond to the power of
% the exponential, for the trigonometric basis the even ones correspond to
% sin and the odd ones correspond to the cosine one.

% parameters
NBases = 11;
NSamples = 1001;

% generate multi-idx to make clear what each element in row is
multi_idx = (-NBases:NBases)';
base_freqs = (-NBases ):1/2:( NBases );

A = zeros( length( multi_idx ), length( base_freqs ) );

for k = 1:length( base_freqs )
    for j = 1:length( multi_idx )
        cidx = base_freqs( k );
        tidx = multi_idx( j )/2;
        if ( mod( multi_idx( j ), 2 ) && ( cidx == tidx || cidx == -tidx ) )
            A( j, k ) = pi;
        elseif(  mod( multi_idx( j ), 2 ) )
            A( j, k ) = ( 2 * tidx * sin( pi * tidx ) * cos( pi * cidx ) )/( tidx^2 - cidx^2 );
        elseif( cidx == tidx )
            A( j, k ) = 1i * pi * sign( cidx * tidx );
        elseif( cidx == -tidx )
            A( j, k ) = 1i * pi * sign( cidx * tidx );
        else
            A( j, k ) = ( -2i * tidx * cos( pi * tidx ) * sin( pi * cidx ) )/( tidx^2 - cidx^2 );
        end
    end
end

% tic
x = linspace( -pi, pi , NSamples );
outF = zeros( NSamples, 1 );
outC = zeros( NSamples, 1 );
% h = A * f1vec';
h = A * boo;
for n = 1:NSamples
    outF( n ) = ( 2 * pi )^-1 * sum( boo .* exp( 1i * x( n ) * base_freqs( : ) ) );
    for m = 1:length( multi_idx )
        if ( mod( multi_idx( m ), 2 ) )
            outC( n ) = outC( n ) + h( m ) * cos( multi_idx( m ) * x( n )/2 );
        else
            outC( n ) = outC( n ) + h( m ) * sin( multi_idx( m ) * x( n )/2 );
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
if (  mod( idx, 2 ) )
    figure; plot( x, cos( x ) .* cos( idx * x/2  ) - idx/2 * sin( x ) .* sin( idx * x/ 2 ) ); hold on;
    %     figure; plot( x, 2 * ( x - 1 ) .* cos( idx * x/2  ) - idx/2 * ( x - 1 ).^2 .* sin( idx * x/ 2 ) ); hold on;
    plot( x, real( outF ), 'r', 'LineWidth', 3 );
    plot( x, real( outC ), 'g' );
else
    figure; plot( x, cos( x ) .* sin( idx * x/2  ) + idx/2 * sin( x ) .* cos( idx * x/ 2 ) ); hold on;
    %     figure; plot( x, 2 * ( x - 1 ) .* sin( idx * x/2  ) + idx/2 * ( x - 1 ).^2 .* cos( idx * x/ 2 ) ); hold on;
    plot( x, real( outF ), 'r', 'LineWidth', 3 );
    plot( x, real( outC ), 'g' );
end

