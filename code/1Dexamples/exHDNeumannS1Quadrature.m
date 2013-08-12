%% parameters
NBases = 102;
compute_nulls = 0;
NSamplesInt = 1501;
NSamplesPlot = 3001;

%% generate multi-idx to make clear what each element in row is
multi_idx = (0:NBases)';

%% compute functions for later use
syms x1
tic
[ t, Weights ] = lgwt( NSamplesInt, -pi, pi );
% defining vector field and its derivative
h1 = matlabFunction( ( x1 - 2 )^4 * ( x1 + 0.5 ) * ( pi^2 - x1^2 ) );
d1 = matlabFunction( simplify( diff( ( x1 - 2 )^4 * ( x1 + 0.5 ) * ( pi^2 - x1^2 ), x1 ) ) );
% in order to get half integer values we scaled the time domain
f1 = h1( t );
f2 = d1( t );
Trig = zeros( length( multi_idx ), NSamplesInt );
DfTrig = zeros( length( multi_idx ), NSamplesInt );
sf = ones( size( multi_idx, 1 ), 1 );
sf( 1 ) = 1/pi;
for k = 1:length( multi_idx )
    if( mod( multi_idx( k ), 2 ) )
        Trig( k, : ) = sin( multi_idx( k ) * t/2 ) * sf( multi_idx( k ) + 1 );
        DfTrig( k, : ) = f1 .* multi_idx( k )/2 .* cos( multi_idx( k ) * t/2 ) + 1/2 * f2 .* sin( multi_idx( k ) * t/2 );
    else
        Trig( k, : ) = cos( multi_idx( k ) * t/2 ) * sf( multi_idx( k ) + 1 );
        DfTrig( k, : ) = f1 .* -multi_idx( k )/2 .* sin( multi_idx( k ) * t/2 ) + 1/2 * f2 .* cos( multi_idx( k ) * t/2 );
    end
end
toc
%% compute the eigenrepresentation of X
% in the complex exponential basis the indices correspond to the power of
% the exponential, for the trigonometric bas
%     figure; plot( 2* atanh( x/pi ), real( squeeze( out( m, : ) ) ) );
%     figure; plot( x, real( squeeze( out( m, : ) ) ) );is the even ones correspond to
% cos and the odd ones correspond to the sin one.
tic
X = zeros( length( multi_idx ) );
parfor k = 1:size( multi_idx, 1 )
    wales = DfTrig( k, : );
    henry = zeros( 1, size( multi_idx, 1 ) );
    for l = 1:size( multi_idx, 1 )
        henry( l ) = sum( Weights' .* wales .* Trig( l, : ) );
    end
    X( k, : ) = henry;
end
Y = X;
for j = 1:size( multi_idx, 1 )
    for k = 1:size( multi_idx, 1 )
        if( multi_idx( j, 1 ) <= NBases - 1 || multi_idx( k, 1 ) <= NBases - 1 )
           Y( j, k ) = 0;
        end
    end
end
toc
%% compute eigenvalues of different portions
if ( compute_nulls )
    tic
    [ ~, helper1 ] = spspaces( X, 2 );
    toc
    V = helper1{1};
    zeroeig = helper1{3}
else
    %         tic
    %         [ V, D ] = eig( X );
    %         toc
    %         eigens = diag(D);
    %         zeroeig = find( abs( ( eigens ) ) <= 1e-4 )
    %         tic
    %             [ U, D, V ] = svd( X );
    %             zeroeig = find( abs( real( diag( D ) ) ) < 1e-10 )
    %             toc
    tic
    zeroeig = 1:2;
    [ V, D ] = eigs( X, length( zeroeig ), 'SM' );
    eigens = diag( D )
    toc;
end
%% just for plotting all the elements inside of the null space
mod_multi_idx = multi_idx;
x = t;
x = linspace( -pi, pi, NSamplesPlot );
% x = linspace( -3, 3, NSamplesPlot );
out = zeros( length( zeroeig ), NSamplesPlot );
parfor m = 1:length( zeroeig )
    for n = 1:NSamplesPlot
        for k = 1:length( mod_multi_idx )
            if ( mod_multi_idx( k ) == 0 ) %% constant term is different (look at the wikipedia on fourier series)
                out( m, n ) = out( m, n ) + V( k, zeroeig( m ) )/sqrt(2);
            elseif ( mod( mod_multi_idx( k ), 2 ) )
                out( m, n ) = out( m, n ) + V( k, zeroeig( m ) ) * sin( mod_multi_idx( k ) * x( n )/2 );
            else
                out( m, n ) = out( m, n ) + V( k, zeroeig( m ) ) * cos( mod_multi_idx( k ) * x( n )/2 );
            end
        end
        out( m, n ) = out( m, n );
    end
end


for m = 1:length( zeroeig )
    blah = squeeze( out( m, : ) ) .* conj( squeeze( out( m, : ) ) );
    figure; plot( x, blah );
    % should be equal to one since the output of V is a normalized vector
    sum( blah ) * ( x( 2 ) - x( 1 ) )/pi
%     sum( Weights' .* blah )/pi
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
