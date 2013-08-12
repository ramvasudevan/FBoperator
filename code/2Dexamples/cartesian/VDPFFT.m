%% parameters
NBases = 1;
compute_nulls = 0;

%% generate multi-idx to make clear what each element in row is
[ gridx, gridy ] = meshgrid( -NBases:NBases );
multi_idx = zeros( size( gridx, 1 ) * size( gridx, 2 ), 2 );
multi_idx( :, 1 ) = gridx( : );
multi_idx( :, 2 ) = gridy( : );

%% compute f_{1,2,3,4} and its DFT
% constructing the signals
Fs = 512;
t = ( -pi ):( ( 2 * pi )/Fs ):( ( pi ) - ( 2 * pi )/Fs );
f1 = zeros( length( t ) ); % x2
f2 = zeros( length( t ) ); % (1 - x1^2)
f3 = zeros( length( t ) ); % (1 - x1^2)*x2 - x1
for j = 1:length( t )
    for k = 1:length( t )
        f1( k, j ) = t( k );
%         f1( k, j ) = 1;
        f2( k, j ) = 1 - t( j )^2;
        f3( k, j ) = ( 1 - t( j )^2 ) * t( k ) - t( j );
    end
end
% computing their FFTs
nfft = length( t );
F1 = fft2( f1, nfft, nfft )/nfft^2;
F2 = fft2( f2, nfft, nfft )/nfft^2;
F3 = fft2( f3, nfft, nfft )/nfft^2;
% phase shifting the FFT
for j = 1:length( t )
    for k = 1:length( t )
        F1( k, j ) = F1( k, j )/( exp( -1i * 2 * pi * Fs/2 * ( k - 1 )/nfft ) ...
            * exp( -1i * 2 * pi * Fs/2 * ( j - 1 )/nfft ) );
        F2( k, j ) = F2( k, j )/( exp( -1i * 2 * pi * Fs/2 * ( k - 1 )/nfft ) ...
            * exp( -1i * 2 * pi * Fs/2 * ( j - 1 )/nfft ) );
        F3( k, j ) = F3( k, j )/( exp( -1i * 2 * pi * Fs/2 * ( k - 1 )/nfft ) ...
            * exp( -1i * 2 * pi * Fs/2 * ( j - 1 )/nfft ) );
    end
end
F1 = fftshift( F1 );
F2 = fftshift( F2 );
F3 = fftshift( F3 );
% frequency vector
freqs = ( -nfft/2:( nfft/2 - 1 ) ) * Fs/nfft;
sampled_freqs = -NBases:NBases;
[ freqsx, freqsy ] = meshgrid( freqs );
f1vec = interp2( freqsx, freqsy, F1, gridx, gridy, 'cubic' );
f2vec = interp2( freqsx, freqsy, F2, gridx, gridy, 'cubic' );
f3vec = interp2( freqsx, freqsy, F3, gridx, gridy, 'cubic' );

%% compute the eigenrepresentation of X
tic
X = zeros( size( multi_idx, 1 ) );
for k = 1:size( multi_idx, 1 )
    c_idx = multi_idx( k, : );
    s_idx( :, 1 ) = c_idx( 1 ) + gridx( : );
    s_idx( :, 2 ) = c_idx( 2 ) + gridy( : );
    f_idx = intersect( intersect( find( c_idx( 1 ) + gridx( : ) >= -NBases ), find( c_idx( 1 ) + gridx( : ) <= NBases ) ), ...
        intersect( find( c_idx( 2 ) + gridy( : ) >= -NBases ), find( c_idx( 2 ) + gridy( : )  <= NBases ) ) );
    h_idx = any( all( bsxfun( @eq, reshape( s_idx.', 1, 2, [] ), multi_idx), 2 ), 3 );
    X( h_idx, k ) = -f2vec( f_idx ) - 1i * ( c_idx( 1 ) * f1vec( f_idx ) + c_idx( 2 ) * f3vec( f_idx ) );
end
toc
%% compute eigenvalues of different portions
if ( compute_nulls )
    [ ~, helper ] = spspaces( X, 2 );
    V = helper{1};
    zeroeig = helper{3};
else
    [ V, D ] = eig( X, eye( size( X ) ) );
    
    zeroeig = [];
    for k = 1:length( D )
        if ( D( k, k ) == 0 )
            zeroeig = [ zeroeig; k ];
        end
    end
end

%% just for plotting
NSamples = 101;
tic
x = linspace( -pi, pi, NSamples );
[ u, v ] = meshgrid( x );
helper = reshape( exp( 1i * kron( u( : ), multi_idx( :, 1 ) ) ).* exp( 1i * kron( v( : ), multi_idx( :, 2 ) ) ), ...
    size( multi_idx, 1 ), NSamples, NSamples );
out = zeros( length( zeroeig ), NSamples, NSamples );
for m = 1:length( zeroeig )
    for n = 1:NSamples
        for l = 1:NSamples
            out( m, n, l ) = V( :, zeroeig( m ) )' * helper( :, l, n );
        end
    end
end
toc
%% actual plotting
for m = 1:size( zeroeig, 2 )
% for m = 75
    figure;
    mesh( x, x, real( squeeze( out( m, :, : ) ) ) ); colorbar;
    pause;
end


