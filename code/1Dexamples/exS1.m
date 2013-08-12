%% parameters
NBases = 101;
compute_nulls = 0;

%% generate multi-idx to make clear what each element in row is
multi_idx = (-NBases:NBases)';

%% compute f_{1-4} and its DFT
syms x1
tic
Fs = 32768; % increase this number to increase the condition number of the matrix (i.e. increase the accuracy of our computation)
%  t = ( -pi ):( ( 2 * pi + 1/Fs )/Fs ):( pi );
t = ( -pi ):( ( 2 * pi )/Fs ):( pi - ( 2 * pi)/Fs );
% t = ( -pi + 5e-10 ):( ( 2 * pi - 1e-9 )/Fs ):( ( pi - 5e-10 ) - ( 2 * pi - 1e-9 )/Fs );
% t = logspace( -pi, pi, Fs );
nfft = length( t );
% defining vector field and its derivative
% h1 = matlabFunction( ( pi^2 - x1^2 )/( 2 * pi ) * -( ( 2 * atanh( x1/pi ) )^2 ) );
% d1 = matlabFunction( -4 * pi * x1/( pi^2 - x1^2 )^2 * ( ( pi^2 - x1^2 )/( 2 * pi ) )^2 * -( ( 2 * atanh( x1/pi ) )^2 )  - ( 2 * ( 2 * atanh( x1/pi ) ) ) );
% f1 = -h1( t );
% f2 = -d1( t );
% f1 = sampvf( t );
% f2 = deriv_sampvf( t );
% f1 = -sampvfS1( t );
% f2 = -deriv_sampvfS1( t );
% h1 = matlabFunction( (pi^2-x1^2)^10/pi^22 * (x1 + 1)^10 * x1 * ( x1 - 2 ) );
% d1 = matlabFunction( diff( (pi^2-x1^2)^10/pi^22 * (x1 + 1)^10 * x1 * ( x1 - 2 ), x1 ) );
% h1 = matlabFunction( (pi^2-x1^2) * ( x1 - 1 ) * ( x1 + 2 ) );
% d1 = matlabFunction( diff( (pi^2-x1^2) * ( x1 - 1 ) * ( x1 + 2 ), x1 ) );
% h1 = matlabFunction( sin( x1 ) );
% d1 = matlabFunction( diff( sin( x1 ), x1 ) );
h1 = matlabFunction( ( x1 - 2 )^2 * ( x1 + 0.5 )^2 * ( pi^2 - x1^2 ) );
d1 = matlabFunction( simplify( diff( ( x1 - 2 )^2 * ( x1 + 0.5 )^2 * ( pi^2 - x1^2 ), x1 ) ) );
f1 = -h1( t );
f2 = -d1( t );
F1 = fftshift( fft( fftshift( f1 ), nfft ) )/nfft * ( 2 * pi );
F2 = fftshift( fft( fftshift( f2 ), nfft ) )/nfft * ( 2 * pi );
% generating actual frequency vector in our preferred bases
freqs = ( -nfft/2:( nfft/2 - 1 ) ) * Fs/nfft;
sampled_freqs = -NBases:NBases;
f1vec = interp1( freqs, F1, sampled_freqs, 'spline' );
f2vec = interp1( freqs, F2, sampled_freqs, 'spline' );
toc
%% compute the eigenrepresentation of X
% 1st coordinate of representation
tic
X = zeros( ( 2 * NBases + 1 ), ( 2 * NBases + 1 )  );
for k = 1:size( multi_idx, 1 )
    c_idx = multi_idx( k );
    s_idx = c_idx + ( -NBases:NBases )';
    h_idx = intersect( find( s_idx >= -NBases ), find( s_idx <= NBases ) );
    
    % actual vector field component of solution
    X( c_idx + h_idx, k ) = ( f1vec( h_idx ) ) * 1i * c_idx + ( f2vec( h_idx ) );
end
toc
%% zero out small entries
% X( abs( X ) < 1e-10 ) = 0;

%% compute eigenvalues of different portions
if ( compute_nulls )
    tic
    [ ~, helper1 ] = spspaces( X, 2 );
    toc
    V = helper1{1};
    zeroeig = helper1{3}
else
    tic
    [ V, D ] = eig( X );
    toc
    eigens = diag(D);
    zeroeig = find( abs( ( eigens ) ) <= 1e-1 )
    %     tic
    %     [ U, D, V ] = svd( X );
    %     zeroeig = find( abs( real( diag( D ) ) ) < 20 )
    %     toc
end

%% just for plotting all the elements inside of the null space
NSamples = 5001;
x = linspace( -pi, pi, NSamples );
out = zeros( length( zeroeig ), NSamples );
for m = 1:length( zeroeig )
    for n = 1:NSamples
        out( m, n ) = (2*pi)^-1 * sum( V( :, zeroeig( m ) ) .* exp( 1i * multi_idx( : ) * x( n ) ) );
    end
end

for m = 1:length( zeroeig )
%     figure; plot( 2* atanh( x/pi ), real( squeeze( out( m, : ) ) ) );
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
