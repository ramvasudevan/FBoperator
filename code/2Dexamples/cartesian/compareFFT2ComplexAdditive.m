B = zeros( ( 2 * NBases + 1 )^2, ( 2 * NBases + 1 )^2 );
helper = zeros( ( 2 * NBases + 1 )^2, 1 );
helper( [ 1 2 ( 2 * NBases + 1 ) + 1 ( 2 * NBases + 1 ) + 2 ] ) = 1;
for k = 1:length( multi_idx )
    B( k:length( multi_idx ), k ) = helper( 1:( length(k:length( multi_idx ) ) ) );
end
A = inv( B );

% tic
x = linspace( -pi, pi , NSamples );
[ u, v ] = meshgrid( x );
[ blahx, blahy ] = meshgrid( NBases:-1:0 );
outF = zeros( NSamples, NSamples );
outS = zeros( NSamples, NSamples );
outC = zeros( NSamples, NSamples );
h = fvec( :, :, 1 );
h = A * h( : );
for n = 1:NSamples
    for l = 1:NSamples
        outF( n, l ) = (4*pi^2)^-1 * sum( sum(  fvec( :, :, 1 ) .* ( exp( 1i * u( n, l ) * gridx ) .* exp( 1i * v( n, l ) * gridy ) ) ) );
        boox = sin( u( n, l ) * blahx );
        boox( blahx == 0 ) = 1;
        booy = sin( v( n, l ) * blahy );
        booy( blahy == 0 ) = 1;
%         outS( n, l ) = (4*pi^2)^-1 * sum( sum(  hvec( :, :, 1 ) .* ( sin( u( n, l ) * blahx ) .* sin( v( n, l ) * blahy ) ) ) );
%         outS( n, l ) = (4*pi^2)^-1 * sum( sum(  hvec( :, :, 1 ) .* ( boox .* booy ) ) );
        outS( n, l ) = -(2048 * (-15 + pi^2) * (-12+pi^2)/pi * sin( v( n, l ) ) * cos( u( n, l )/2 ) );
        for m = 1:length( h )
            outC( n, l ) = outC( n, l ) + h( m ) * sum( exp( 1i * u( n, l ) * multi_idx( find( B( :, m  ) ), 1  ) ) .* exp( 1i * v( n, l ) * multi_idx( find( B( :, m  ) ), 2  ) ) );
        end
        outC( n, l ) = outC( n, l ) * 0.25 * (4*pi^2)^-1;
    end
end
% toc

