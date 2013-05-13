%% parameters
NBases = 10;
compute_nulls = 0;

%% generate multi-idx to make clear what each element in row is
multi_idx = zeros( ( 2 * NBases + 1 )^2, 2 );
counter = 0;
for m = -NBases:NBases
    for n = -NBases:NBases
        counter = counter + 1;
        multi_idx( counter, : ) = [ m n ];
    end
end

%% compute f_{1,2,3,4} and its DFT
% constructing the signals
Fs = 1000;
t = ( -pi ):( ( 2 * pi )/Fs ):( pi );
f1 = -t;
f2 = t;
f3 = ( t.^2 - 1 );
f4 = zeros( length( t ) );
for j = 1:length( t )
    for k = 1:length( t )
        f4( j, k ) = ( t( j )^2 - 1 ) * t( k );
    end
end
% computing their FFTs
nfft = length( t );
F1 = fftshift( fft( f1, nfft ) );
F2 = fftshift( fft( f2, nfft ) );
F3 = fftshift( fft( f3, nfft ) );
F4 = fftshift( fft2( f4, nfft, nfft ) );
% frequency vector
freqs = ( -nfft/2:( nfft/2 - 1 ) ) * Fs/nfft;
sampled_freqs = -NBases:NBases;
f1vec = interp1( freqs, F1, sampled_freqs );
f2vec = interp1( freqs, F2, sampled_freqs );
f3vec = interp1( freqs, F3, sampled_freqs );
[ freqsx, freqsy ] = meshgrid( freqs );
[ sampled_freqsx, sampled_freqsy ] = meshgrid( sampled_freqs );
f4vec = interp2( freqsx, freqsy, F4, sampled_freqsx, sampled_freqsy );
%% compute the eigenrepresentation of X
% 1st coordinate of representation
tic
X1 = zeros( ( 2 * NBases + 1 )^2, ( 2 * NBases + 1 )^2  );
for k = 1:size( multi_idx, 1 )
    c_idx = multi_idx( k, : );
    s_idx = zeros( 2 * NBases + 1, 2 );
    s_idx( :, 1 ) = c_idx( 1 );
    s_idx( :, 2 ) = c_idx( 2 ) + ( -NBases:NBases )';
    h_idx = any( all( bsxfun( @eq, reshape( s_idx.', 1, 2, [] ), multi_idx), 2 ), 3 );
    X1( h_idx, k ) = -1i * c_idx( 1 ) * f1vec( intersect( find( c_idx( 2 ) + (-NBases:NBases) >= -NBases ), find( c_idx( 2 ) + (-NBases:NBases) <= NBases ) ) );
%     for j = -NBases:NBases
%         h_idx = ismember( multi_idx, [ c_idx( 1 ) c_idx( 2 ) + j ], 'rows' );
% %         h_idx = any( all( bsxfun( @eq, reshape( [ c_idx( 1 ) c_idx( 2 ) + j ].', 1, 2, [] ), multi_idx), 2 ), 3 );
%         X1( h_idx, k ) = -1i * c_idx( 1 ) * f1vec( j + NBases + 1 );
%     end
end
toc

% 2nd coordinate of representation
tic
X2 = zeros( ( 2 * NBases + 1 )^2, ( 2 * NBases + 1 )^2  );
for k = 1:size( multi_idx, 1 )
    c_idx = multi_idx( k, : );
    for j = -NBases:NBases
        h_idx = ismember( multi_idx, [ c_idx( 1 ) + j c_idx( 2 ) ], 'rows' );
        X2( h_idx, k ) = 1i * c_idx( 2 ) * f2vec( j + NBases + 1 )  + f3vec( j + NBases + 1 );
        
        for l = -NBases:NBases
            h_idx = ismember( multi_idx, [ c_idx( 1 ) + j c_idx( 2 ) + l ], 'rows' );
            X2( h_idx, k ) = X2( h_idx, k ) + 1i * c_idx( 2 ) * f4vec( j + NBases + 1, l + NBases + 1 );
        end
    end
end
toc
%% compute eigenvalues of different portions
if ( compute_nulls )
    [ ~, helper1 ] = spspaces( X1, 2 );
    V1 = helper1{1};
    zeroeig1 = helper1{3};
    [ ~, helper2 ] = spspaces( X2, 2 );
    V2 = helper2{1};
    zeroeig2 = helper2{3};
else
    [ V1, D1 ] = eig( X1, eye( size( X1 ) ) );
    [ V2, D2 ] = eig( X2, eye( size( X2 ) ) );
    
    zeroeig1 = [];
    zeroeig2 = [];
    for k = 1:length( D1 )
        if ( D1( k, k ) == 0 )
            zeroeig1 = [ zeroeig1; k ];
        end
        if ( D2( k, k ) == 0 )
            zeroeig2 = [ zeroeig2; k ];
        end
    end
end


%% just for plotting
NSamples = 101;
x = linspace( -pi, pi, NSamples );
outt = zeros( length( zeroeig1 ), NSamples, NSamples );
for m = 1:length( zeroeig1 )
    for n = 1:NSamples
        for j = 1:NSamples
            outt( m, n, j ) = outt( m, n, j ) + V1( :, zeroeig1( m ) )' * exp( 1i * multi_idx( :, 1 ) * x( n ) + 1i *  multi_idx( :, 2 ) * x( j ) );
        end
    end
end

figure; hold on;
for m = 1:length( zeroeig1 )
    mesh( x, x, real( squeeze( outt( m, :, : ) ) ) );
end


