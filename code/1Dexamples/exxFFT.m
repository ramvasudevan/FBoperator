%% parameters
NBases = 1;
compute_nulls = 1;

%% generate multi-idx to make clear what each element in row is
multi_idx = (-NBases:NBases)';

%% compute f_{1,2} and its DFT
Fs = 512;
% t = 0:( ( 2 * pi )/Fs ):( ( 2 * pi ) - ( 2 * pi)/Fs );
t = ( -pi ):( ( 2 * pi )/Fs ):( ( pi ) - ( 2 * pi )/Fs );
f1 = ones( size( t ) ); % I
f1(  t <= -pi/2  ) = 0;
f1(  t > pi/2  ) = 0;
f4 = zeros( size( t ) ); % Ix
f4( intersect( find( t >= -pi/2 ), find( t <= pi/2 ) ) ) = t( intersect( find( t >= -pi/2 ), find( t <= pi/2 ) ) );
% computing their FFTs
nfft = length( t );
freqs = ( -nfft/2:( nfft/2 - 1 ) ) * Fs/nfft;
F1 = fft( f1, nfft )/nfft * ( 2 * pi );
F2 = -pi/2 * exp( -1i * pi/2 * freqs ); % delta function * x at -pi/2
F3 = -pi/2 * exp( 1i * pi/2 * freqs ); % detla function * x at pi/2
F4 = fft( f4, nfft )/nfft * ( 2 * pi );
% phase shifting the FFT
% for j = 1:length( t )
%     F1( j ) = F1( j )/( exp( -1i * 2 * pi * Fs/2 * ( j - 1 )/nfft ) );
%     F2( j ) = F2( j )/( exp( -1i * 2 * pi * Fs/2 * ( j - 1 )/nfft ) );
% end
F1 = fftshift( F1 );
F4 = fftshift( F4 );
% frequency vector for analysis
sampled_freqs = -NBases:NBases;
f1vec = interp1( freqs, F1, sampled_freqs, 'cubic' );
f2vec = interp1( freqs, F2, sampled_freqs, 'cubic' );
f3vec = interp1( freqs, F3, sampled_freqs, 'cubic' );
f4vec = interp1( freqs, F4, sampled_freqs, 'cubic' );
%% compute the eigenrepresentation of X
% 1st coordinate of representation
tic
X= zeros( ( 2 * NBases + 1 ), ( 2 * NBases + 1 )  );
for k = 1:size( multi_idx, 1 )
    c_idx = multi_idx( k );
    %     for j = -NBases:NBases
    %         if ( c_idx + j >= -NBases && c_idx + j <= NBases )
    %             X1( c_idx + j + NBases + 1, k ) = f1vec( j + NBases + 1 ) + f2vec( j + NBases + 1 ) * 1i * c_idx( 1 );
    %         end
    %     end
    s_idx = c_idx + ( -NBases:NBases )';
    h_idx = intersect( find( s_idx >= -NBases ), find( s_idx <= NBases ) );
    X( c_idx + h_idx, k ) = ( f1vec( h_idx )+ f2vec( h_idx ) + f3vec( h_idx ) + ...
        f4vec( h_idx ) * 1i * c_idx( 1 ) );
end
toc
%% zero out small entries
% X1( abs( X1 ) < 1e-10 ) = 0;

%% compute eigenvalues of different portions
if ( compute_nulls )
    tic
    [ ~, helper1 ] = spspaces( X, 2 );
    toc
    V = helper1{1};
    zeroeig = helper1{3};
else
    tic
    [ V, D ] = eig( X, eye( size( X ) ) );
    toc
    zeroeig = [];
    for k = 1:length( D )
        if ( abs( D( k, k ) ) == 0 )
            zeroeig = [ zeroeig; k ];
        end
    end
end

% %% just for plotting all the elements inside of the null space
% NSamples = 5001;
% x = linspace( -pi, pi, NSamples );
% out = zeros( length( zeroeig ), NSamples );
% for m = 1:length( zeroeig )
%     for n = 1:NSamples
%         out( m, n ) = out( m, n ) + V( :, zeroeig( m ) )' * exp( 1i * multi_idx( : ) * x( n ) );
%     end
% end
% 
% figure; hold on;
% for m = 1:length( zeroeig )
%     plot( x, real( squeeze( out( m, : ) ) ) );
% end

