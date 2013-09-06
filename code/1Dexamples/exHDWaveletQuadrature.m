%% parameters
WaveletScale = 1;
NBases = 2^(WaveletScale);
compute_nulls = 0;
NSamplesWavelet = 8; % its actually going to have resolution 2^(-NSamplesWavelet)
NSamplesInt = 15001;
addpath('ex_vf');

%% generate base scaling function and its derivative
scalingFilter = daub( 2 * 10 );
[ tWavelet, phi, ~ ] = phivals( scalingFilter, NSamplesWavelet, 0 );
[ tWavelet, dphi, ~ ] = phivals( 2 * scalingFilter, NSamplesWavelet, 1 );

%% generate offset to make clear what each element in row is
offset_max = floor( 2^( WaveletScale ) * 5 );
offset_min = ceil( -2^( WaveletScale ) * 5 - tWavelet( end ) );
multi_idx = ( offset_min:offset_max )';

%% generate quadrature formulas
[ t, Weights ] = lgwt( NSamplesInt, -5, 5 );
% t = -5:( 2^( -NSamplesWavelet ) * 2^( -WaveletScale ) ):5;
% Weights = t( 2 ) - t( 1 );

%% compute functions for later use
syms x1
tic
% defining vector field and its derivative
h1 = matlabFunction( ( x1 + 3 ) * ( x1 - 3 )^2 * ( 5^2 - x1^2 ) );
d1 = matlabFunction( simplify( diff( ( x1 + 3 ) * ( x1 - 3 )^2 * ( 5^2 - x1^2 ), x1 ) ) );
% precompute several quantities that are repeatedly used
f1 = h1( t );
f2 = d1( t );
% renormalize wavelet to appropriate scale and [ 0, 2^( -WaveletScale ) * tWavelet(end) ]
tMRA = tWavelet * 2^( -WaveletScale );
PhiBases = 2^( WaveletScale/2 ) * phi;
dPhiBases = 2^( WaveletScale/2 ) * 2^( WaveletScale ) * dphi;
% fill in arrays
% Trig = zeros( length( multi_idx ), NSamplesInt );
% DfTrig = zeros( length( multi_idx ), NSamplesInt );
Trig = zeros( length( multi_idx ), length( t ) );
DfTrig = zeros( length( multi_idx ), length( t ) );
for k = 1:length( multi_idx )
    tBases = ( tWavelet + multi_idx( k ) ) * 2^( -WaveletScale );
    %     tBases = ( tWavelet * 2^( -WaveletScale ) - 2^( -WaveletScale ) * multi_idx( k ) );
    Trig( k, : ) = interp1( tBases, PhiBases, t, 'linear', 'extrap' );
    helperD = interp1( tBases, dPhiBases, t, 'linear', 'extrap' );
    Trig( k, t <= tBases( 1 ) ) = 0;
    Trig( k, t >= tBases( end ) ) = 0;
    helperD( t <= tBases( 1 ) ) = 0;
    helperD( t >= tBases( end ) ) = 0;
%     Trig( k, isnan( Trig( k, : ) ) ) = 0;
%     helperD( isnan( helperD ) ) = 0;
    DfTrig( k, : ) = f1 .* helperD + 1/2 .* f2 .* Trig( k, : )';
end
toc
%% compute the eigenrepresentation of X
% in the complex exponential basis the indices correspond to the power of
% the exponential, for the trigonometric basis the even ones correspond to
% cos and the odd ones correspond to the sin one.
tic
X = zeros( length( multi_idx ) );
for k = 1:size( multi_idx, 1 )
    wales = DfTrig( k, : );
    henry = zeros( 1, size( multi_idx, 1 ) );
    for l = 1:size( multi_idx, 1 )
        henry( l ) = sum( Weights' .* wales .* Trig( l, : ) );
    end
    X( k, : ) = henry;
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
    tic
    [ V, D ] = eig( X );
    toc
    eigens = diag(D);
    zeroeig = find( abs( ( eigens ) ) <= 1e-4 )
    %         tic
    %             [ U, D, V ] = svd( X );
    %             zeroeig = find( abs( real( diag( D ) ) ) < 1e-10 )
    %             toc
    %     tic
    %     zeroeig = 1:2;
    %     [ V, D ] = eigs( X, length( zeroeig ), 'SM' );
    %     eigens = diag( D )
    %     toc;
end
%% just for plotting all the elements inside of the null space
x = t;
% out = V( :, zeroeig )' * Trig;
% for m = 1:length( zeroeig )
%     blah = squeeze( out( m, : ) ) .* conj( squeeze( out( m, : ) ) );
%     figure; plot( x, blah );
% end

%% plot vector field
[ T, Z ] = ode45( @polyR1_plot, [ 0, 1 ], 0  );
figure; plot( T, Z );

%% plot mesa function evolution
syms x1 w1;
T = 0.05;
h8 = @(x1) 10 * (heaviside( x1 + 1.2 ) - heaviside( x1 + 0.8 ) );
% f8 = h8( t );
f8 = smooth( h8( t ), 1000, 'loess' );
% determine the original coefficients
mesaout0 = zeros( size( multi_idx ) );
for i = 1:size( Trig, 1 )
    mesaout0( i ) = sum( Weights .* f8 .* Trig( i, : )' );
end
% mesaout0( abs( mesaout0 ) < 1e-3 ) = 0;
% X( abs( X ) < 1000 ) = 0;

% maybe worth looking at the expm_conditioning of the result to ensure
% that it looks okay, also there are faster and stabler ways to compute
% this (balance and then exponentiate...
stupid = linspace( 0, T, 10 );
for k = 1:length( stupid )
    expmXp = expm_new( X * stupid( k ) );
%     expmXp( abs( expmXp ) < 1e-1 ) = 0;
    figure; spy( expmXp );
%     expmXp = V * expm( D * stupid( k ) ) * V';
    helper8 = expmXp * mesaout0;
    mesaout8 = helper8' * Trig;
    h = figure; hold on;
    plot( x, mesaout8 .* conj( mesaout8 ) );
    title( sprintf( 't = %f', stupid( k ) ) );
    %     filename = sprintf( 'indicator_mesa_images/image_%2.3f.pdf', stupid( k ) );
    %     print( h, '-dpdf', filename );
    pause;
end