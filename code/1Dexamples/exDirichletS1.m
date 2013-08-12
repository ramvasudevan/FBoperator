%% parameters
NBases = 11;
compute_nulls = 0;
NSamples = 1001;

%% generate multi-idx to make clear what each element in row is
multi_idx = (-NBases:NBases)';

%% compute f_{1-4} and its DFT
syms x1
tic
Fs = 16384; % increase this number to increase the condition number of the matrix (i.e. increase the accuracy of our computation)
t = ( -pi ):( ( 2 * pi )/Fs ):( pi - ( 2 * pi)/Fs );
nfft = length( t );
% defining vector field and its derivative
% h1 = matlabFunction( ( x1 - 1 )^2 );
% d1 = matlabFunction( diff( ( x1 - 1 )^2, x1 ) );
h1 = matlabFunction( sin( x1 ) );
d1 = matlabFunction( diff( sin( x1 ), x1 ) );
% in order to get half integer values we scaled the time domain
f1 = h1( 2 * t );
f2 = d1( 2 * t );
F1 = fftshift( fft( fftshift( f1 ), nfft ) )/nfft * ( 2 * pi );
F2 = fftshift( fft( fftshift( f2 ), nfft ) )/nfft * ( 2 * pi );
% generating actual frequency vector in our preferred bases
freqs = ( -nfft/2:( nfft/2 - 1 ) ) * Fs/nfft;
% remember in order to get the half integer values we scaled the signal up
sampled_freqs = (-4*NBases):(4*NBases); 
f1vec = interp1( freqs, F1, sampled_freqs, 'spline' );
f2vec = interp1( freqs, F2, sampled_freqs, 'spline' );
% we relabel the constructed signal back to its half integer values
sampled_freqs = (-2 * NBases ):1/2:( 2 * NBases );
base_freqs = (-NBases ):1/2:( NBases );
% f1vec( mod( sampled_freqs, 1 ) ~= 0  ) = 0;
% f2vec( mod( sampled_freqs, 1 ) ~= 0 ) = 0;
toc
%% Construct Projection onto our Trig bases
% in the complex exponential basis the indices correspond to the power of
% the exponential, for the trigonometric basis the even ones correspond to
% sin and the odd ones correspond to the cosine one.
tic
P = zeros( length( multi_idx ), length( base_freqs ) );
for k = 1:length( base_freqs )
    for j = 1:length( multi_idx )
        cidx = base_freqs( k );
        tidx = multi_idx( j )/2;
        if ( mod( multi_idx( j ), 2 ) && ( cidx == tidx || cidx == -tidx ) )
            P( j, k ) = pi;
        elseif(  mod( multi_idx( j ), 2 ) )
            P( j, k ) = ( 2 * tidx * sin( pi * tidx ) * cos( pi * cidx ) )/( tidx^2 - cidx^2 );
        elseif( cidx == tidx )
            P( j, k ) = 1i * pi * sign( cidx * tidx );
        elseif( cidx == -tidx )
            P( j, k ) = 1i * pi * sign( cidx * tidx );
        else 
            P( j, k ) = ( -2i * tidx * cos( pi * tidx ) * sin( pi * cidx ) )/( tidx^2 - cidx^2 );
        end
    end
end
toc
%% compute the eigenrepresentation of X
% 1st coordinate of representation
tic
X = zeros( ( 2 * NBases + 1 ) );
for k = 1:size( multi_idx, 1 )
    helper = zeros( length( base_freqs ), 2 );
    c_idx = multi_idx( k );
    
    if( mod( c_idx, 2 ) )
        % cos( kx/2 ) k is odd case
        s_idx = c_idx/2 + sampled_freqs';
        h_idx = intersect( find( s_idx >= -NBases ), find( s_idx <= NBases ) );
        helper( :, 1 ) = 1/2 * f2vec( h_idx ) + 1i * c_idx/4 * f1vec( h_idx );
        
        s_idx = -c_idx/2 + sampled_freqs';
        h_idx = intersect( find( s_idx >= -NBases ), find( s_idx <= NBases ) );
        helper( :, 2 ) = 1/2 * f2vec( h_idx ) - 1i * c_idx/4 * f1vec( h_idx );
    else
        % sin( kx/2 ) k is even case
        s_idx = c_idx/2 + sampled_freqs';
        h_idx = intersect( find( s_idx >= -NBases ), find( s_idx <= NBases ) );
        helper( :, 1 ) = -1i * 1/2 * f2vec( h_idx ) + c_idx/4 * f1vec( h_idx );
        
        s_idx = -c_idx/2 + sampled_freqs';
        h_idx = intersect( find( s_idx >= -NBases ), find( s_idx <= NBases ) );
        helper( :, 2 ) = 1i * 1/2 * f2vec( h_idx ) + c_idx/4 * f1vec( h_idx );
    end
    
%     if ( multi_idx( k ) == -1 )
%         boo = ( helper( :, 1 ) + helper( :, 2 ) );
%     end

    if ( k >= 3 && k <= length( multi_idx ) - 3 )
        compareFFT2Trig( multi_idx( k ), helper( :, 1 ) + helper( :, 2 ) );
    end
    
    % actual vector field component of solution
    X( :, k ) = P * ( helper( :, 1 ) + helper( :, 2 ) );
end
toc
X = X( [ 1:NBases ( NBases + 2 ):( 2 * NBases + 1 ) ], : );
X = X( :, [ 1:NBases ( NBases + 2 ):( 2 * NBases + 1 ) ] );
X = X( 3:( length( multi_idx ) - 3 ), 3:( length( multi_idx ) - 3 ) );

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
    zeroeig = find( abs( ( eigens ) ) <= 5e-15 )
%     tic
%         [ U, D, V ] = svd( X );
%         zeroeig = find( abs( real( diag( D ) ) ) < 20 )
%         toc
%         tic
    %     zeroeigR = 1:5;
    %     [ VR, D ] = eigs( X, length( zeroeigR ), 'SM' );
    %     eigens = diag( D );
    %     toc;
end

%% just for plotting all the elements inside of the null space
mod_multi_idx = [ ( -NBases + 2 ):( -1 ) ( 1:( NBases - 2 ) ) ];
% mod_multi_idx = [ ( -NBases):( -1 ) ( 1:( NBases) ) ];
x = linspace( -pi, pi, NSamples );
out = zeros( length( zeroeig ), NSamples );
for m = 1:length( zeroeig )
    for n = 1:NSamples
        for k = 1:length( mod_multi_idx )
            if ( mod( mod_multi_idx( k ), 2 ) )
                out( m, n ) = out( m, n ) + V( k, zeroeig( m ) ) * cos( mod_multi_idx( k ) * x( n )/2 );
            else
                out( m, n ) = out( m, n ) + V( k, zeroeig( m ) ) * sin( mod_multi_idx( k ) * x( n )/2 );
            end
        end
        out( m, n ) = out( m, n ) * ( 2 * pi )^-2;
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
