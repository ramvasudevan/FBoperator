%% parameters
NBases = 101;
compute_nulls = 1;

%% generate multi-idx to make clear what each element in row is
multi_idx = (-NBases:NBases)';

%% compute f_{1-4} and its DFT
syms x1
tic
Fs = 4096; % increase this number to increase the condition number of the matrix (i.e. increase the accuracy of our computation)
t = -pi:( ( 2 * pi )/Fs ):( pi - ( 2 * pi)/Fs );
nfft = length( t );
% defining vector field and its derivative
% h1 = matlabFunction( exp(-(x1/3)^100) * -x1 );
% d1 = matlabFunction( diff( exp(-(x1/3)^100) * -x1, x1 ) );
h1 = matlabFunction( exp( -x1^100 ) * -x1 );
d1 = matlabFunction( diff( exp( -x1^100 ) * -x1, x1 ) );
% h1 = @(x) -x;
% d1 = @(x) -ones(size(x));
% building vectors to FFT
f1 = -h1( t );
f2 = -d1( t );
F1 = fftshift( fft( fftshift( f1 ), nfft ) )/nfft * ( 2 * pi );
F2 = fftshift( fft( fftshift( f2 ), nfft ) )/nfft * ( 2 * pi );
% generating actual frequency vector in our preferred bases
freqs = ( -nfft/2:( nfft/2 - 1 ) ) * Fs/nfft;
sampled_freqs = -NBases:NBases;
f1vec = interp1( freqs, F1, sampled_freqs, 'spline' );
f2vec = interp1( freqs, F2, sampled_freqs, 'spline' );
% f1vec = zeros(size(sampled_freqs));
% f2vec = zeros(size(sampled_freqs));
% FFT of boundary components
f3vec = h1( 1 ) * exp( 1i * -sampled_freqs * 1 )/( 2 * pi );
f4vec = -h1( -1 ) * exp( 1i * sampled_freqs * 1 )/( 2 * pi );
toc
%% compute the eigenrepresentation of X
% 1st coordinate of representation
tic
X= zeros( ( 2 * NBases + 1 ), ( 2 * NBases + 1 )  );
for k = 1:size( multi_idx, 1 )
    c_idx = multi_idx( k );
    s_idx = c_idx + ( -NBases:NBases )';
    h_idx = intersect( find( s_idx >= -NBases ), find( s_idx <= NBases ) );
    
    % actual vector field component of solution
    X( c_idx + h_idx, k ) = ( f1vec( h_idx ) ) * 1i * c_idx + ( f2vec( h_idx ) );
    
    % actual boundary components of solution
    X( :, k ) = X( :, k ) + ( f3vec * exp( 1i * c_idx * 1 ) + f4vec * exp( -1i * c_idx * 1 ) )';
%     X( :, k ) = f3vec + f4vec;
    %X( c_idx + h_idx, k ) = ( f1vec( h_idx ) + f3vec( h_idx ) ) * 1i * c_idx + ( f2vec( h_idx ) + f4vec( h_idx ) );
end
toc
%% zero out small entries
% X( abs( X ) < 1e-10 ) = 0;

%% compute eigenvalues of different portions
if ( compute_nulls )
    tic
    [ ~, helper1 ] = spspaces( X, 2 , 5e-2 );
    toc
    V = helper1{1};
    zeroeig = helper1{3};
else
    tic
    [ V, D ] = eig( X, eye( size( X ) ) );
    toc
    eigens = diag(D);
    zeroeig = find( abs( eigens ) <= 1e-3 );
end

%% just for plotting all the elements inside of the null space
NSamples = 1001;
x = linspace( -pi, pi, NSamples );
out = zeros( length( zeroeig ), NSamples );
for m = 1:length( zeroeig )
    for n = 1:NSamples
        out( m, n ) = (2*pi)^-1 * sum( V( :, zeroeig( m ) ) .* exp( 1i * multi_idx( : ) * x( n ) ) );
    end
end

for m = 1:length( zeroeig )
    figure; plot( x, real( squeeze( out( m, : ) ) ) );
end

%% plot mesa function evolution
syms x1 w1;
T = 1;
h8 = @(x1) 1 * (heaviside( x1 + 0.2 ) - heaviside( x1 - 0.2 ) );
f8 = h8( t );
F8 = fftshift( fft( fftshift( f8 ), nfft ) )/nfft;
xvec8 = interp1( freqs, F8, sampled_freqs, 'nearest' );
xvec8 = reshape( xvec8, size( xvec8, 1 ) * size( xvec8, 2 ), 1 );
% xvec = ones( size( multi_idx( :, 1 ) ) );
% xvec( isnan( xvec ) ) = 0; xvec( isinf( xvec ) ) = 1e10;

% maybe worth looking at the expm_conditioning of the result to ensure
% that it looks okay, also there are faster and stabler ways to compute
% this (balance and then exponentiate...
stupid = linspace( 0, T, 20 );
for k = 1:length( stupid )
    expmXp = expm( X * stupid( k ) );
    helper8 = expmXp * xvec8;
    % NSamples = 201;
    x = linspace( -pi, pi, NSamples );
    mesaout8 = zeros( NSamples, 1 );
    for n = 1:NSamples
        mesaout8( n ) = (2*pi)^-1 * sum( helper8 .* ( exp( 1i * x( n ) * multi_idx( : ) ) ) );
    end
    h = figure; hold on;
    plot( x, real( mesaout8 ) );
%     filename = sprintf( 'indicator_mesa_images/image_%2.3f.pdf', stupid( k ) );
%     print( h, '-dpdf', filename );
    pause;
end
