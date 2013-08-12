%% parameters
NBases = 11;
compute_nulls = 1;
scaling = 3;

%% generate multi-idx to make clear what each element in row is
[ gridx, gridy ] = meshgrid( -NBases:NBases );
multi_idx = zeros( size( gridx, 1 ) * size( gridx, 2 ), 2 );
multi_idx( :, 1 ) = gridx( : );
multi_idx( :, 2 ) = gridy( : );

%% compute f_{1-5} and its DFT
% constructing the signals
syms x1 x2;
tic
Fs = 4096; % increase this number to increase the condition number of the matrix (i.e. increase the accuracy of our computation)
t = -pi:( ( 2 * pi )/Fs ):( pi - ( 2 * pi)/Fs );
% t = 0:( ( 2 * pi )/Fs ):( 2 * pi - ( 2 * pi)/Fs );
[ tx, ty ] = meshgrid( t );
h1 = matlabFunction( scaling * x2 * exp( -(x1/pi)^100 ) * exp( -(x2/pi)^100 ) );
h2 = matlabFunction( -(scaling*x1 + scaling*x2.*((scaling*x1).^2 - 1)) * exp( -(x1/pi)^100 ) * exp( -(x2/pi)^100 ) );
d1 = matlabFunction( diff( scaling * x2 * exp( -(x1/pi)^100 ) * exp( -(x2/pi)^100 ), x1 ) );
d2 = matlabFunction( diff( -(scaling*x1 + scaling*x2.*((scaling*x1).^2 - 1)) * exp( -(x1/pi)^100 ) * exp( -(x2/pi)^100 ), x2 ) );
% h1 = @(x1,x2) scaling * x2 .* exp(-1./(1-(x1/pi).^100)) .* exp(-1./(1-(x2/pi).^100)) * (exp(1))^2;
% h2 = @(x1,x2) -(scaling*x1 + scaling*x2.*((scaling*x1).^2 - 1)) .* exp(-1./(1-(x1/pi).^100)) .* exp(-1./(1-(x2/pi).^100)) * (exp(1))^2;
% d1 = @(x1,x2) scaling * x2 .* -(100*pi^100*exp(-1./(1-x1.^100/pi^100)).*x1.^99)./(pi^100-x1.^100).^2 .* exp(-1./(1-(x2/pi).^100)) * (exp(1))^2;
% d2 = @(x1,x2) -( scaling * exp(-1./(1-x2.^100/pi^100)).*(x2.^200.*(scaling^2*x1.^2-1)-2*pi^100*x2.^99.*(51*x2.*(scaling^2*x1.^2-1)+50*x1)+pi^200.*(scaling^2*x1.^2-1)))./(pi^100-x2.^100).^2 .* exp(-1./(1-(x1/pi).^100)) .* (exp(1))^2;
f1 = -h1( tx, ty );
f2 = -d1( tx, ty );
f3 = -h2( tx, ty );
f4 = -d2( tx, ty );
% f1( isnan( f1 ) ) = 0; f1( isinf( f1 ) ) = 0;
% f2( isnan( f2 ) ) = 0; f2( isinf( f2 ) ) = 0;
% f3( isnan( f3 ) ) = 0; f3( isinf( f3 ) ) = 0;
% f4( isnan( f4 ) ) = 0; f4( isinf( f4 ) ) = 0;
% computing their FFTs
nfft = length( t );
F1 = fftshift( fft2( fftshift( f1 ), nfft, nfft ) )/nfft^2;
F2 = fftshift( fft2( fftshift( f2 ), nfft, nfft ) )/nfft^2;
F3 = fftshift( fft2( fftshift( f3 ), nfft, nfft ) )/nfft^2;
F4 = fftshift( fft2( fftshift( f4 ), nfft, nfft ) )/nfft^2;
% frequency vector
freqs = ( -nfft/2:( nfft/2 - 1 ) ) * Fs/nfft;
sampled_freqs = -NBases:NBases;
[ freqsx, freqsy ] = meshgrid( freqs );
f1vec = interp2( freqsx, freqsy, F1, gridx, gridy, 'cubic' );
f2vec = interp2( freqsx, freqsy, F2, gridx, gridy, 'cubic' );
f3vec = interp2( freqsx, freqsy, F3, gridx, gridy, 'cubic' );
f4vec = interp2( freqsx, freqsy, F4, gridx, gridy, 'cubic' );
toc
%% compute the eigenrepresentation of X
tic
Xnum = zeros( size( multi_idx, 1 ) );
Z = zeros( size( multi_idx, 1 ) );
s_idx = zeros( size( gridx, 1 ) * size( gridx, 2 ), 2 );
for k = 1:size( multi_idx, 1 )
    c_idx = multi_idx( k, : );
    s_idx( :, 1 ) = c_idx( 1 ) + gridx( : );
    s_idx( :, 2 ) = c_idx( 2 ) + gridy( : );
    f_idx = intersect( intersect( find( c_idx( 1 ) + gridx( : ) >= -NBases ), find( c_idx( 1 ) + gridx( : ) <= NBases ) ), ...
        intersect( find( c_idx( 2 ) + gridy( : ) >= -NBases ), find( c_idx( 2 ) + gridy( : )  <= NBases ) ) );
    h_idx = any( all( bsxfun( @eq, reshape( s_idx.', 1, 2, [] ), multi_idx), 2 ), 3 );
%     X( h_idx, k ) = 1i * c_idx( 1 ) * f1vec( f_idx ) + f2vec( f_idx ) + f3vec( f_idx ) + ...
%         1i * c_idx( 2 ) * f4vec( f_idx ) + f5vec( f_idx ) + f6vec( f_idx );
    Xnum( h_idx, k ) = 1i * c_idx( 1 ) * f1vec( f_idx ) + f2vec( f_idx ) + ...
        1i * c_idx( 2 ) * f3vec( f_idx ) + f4vec( f_idx );
end
toc
%% compute eigenvalues of different portions
if ( compute_nulls )
    [ helperL, helperR ] = spspaces( Xnum, 3, 10 );
    VR = helperR{1};
    zeroeigR = helperR{3};
    VL = helperL{1};
    zeroeigL = helperL{3};
else
    [ VR, D ] = eig( Xnum, eye( size( Xnum ) ) );
    
    zeroeigR = find( abs( real( diag( D ) ) ) < 1e-3 );
end

%% just for plotting
tic
NSamples = 101;
x = linspace( -pi, pi, NSamples );
[ u, v ] = meshgrid( x );
z_idx = any( all( bsxfun( @eq, reshape( [ 0 0 ].', 1, 2, [] ), multi_idx), 2 ), 3 );
outR = zeros( length( zeroeigR ), NSamples, NSamples );
% outL = zeros( length( zeroeigL ), NSamples, NSamples );
for m = 1:length( zeroeigR )
    for n = 1:NSamples
        for l = 1:NSamples
            outR( m, n, l ) = (4*pi^2)^-1 * sum( VR( :, zeroeigR( m ) ) .* ( exp( 1i * u( n, l ) * multi_idx( :, 1 ) ) .* exp( 1i * v( n, l ) * multi_idx( :, 2 ) ) ) ) ;        
        end
    end
end
toc
%% actual plotting
for m = 1:length( zeroeigR )
    % for m = 75
    figure;
    surf( u, v, real( squeeze( outR( m, :, : ) ) ) ); colorbar;
    xlabel('x1');
    ylabel('x2');
    shading interp;
    figure;
    contour( u, v, imag( squeeze( outR( m, :, : ) ) ) ); colorbar;
    xlabel('x1');
    ylabel('x2');
    shading interp;
end

%% plot vector field
[ u, v ] = meshgrid( -pi:0.05*pi:( pi ), -pi:0.05*pi:( pi ) );
Du = h1( u, v );
Dv = h2( u, v );
figure; hold on;
quiver( u, v, Du, Dv );
[ ~, Y ] = ode45( @vdpvf, [0 1000], [0 1.2] );
plot( Y( :, 1 ), Y( :, 2 ), 'g' );
% [ ~, Y ] = ode45( @vdpvf, [0 1000], [ 0  pi-1.5e-1 ] );
% plot( Y( :, 1 ), Y( :, 2 ), 'r' );
%% plot mesa function evolution
syms x1 x2 w1 w2;
T = 0;
F6 = matlabFunction( simplify( fourier( fourier( (heaviside(x1+1.5) - heaviside(x1+0.5))*(heaviside(x2+1/2) - heaviside(x2-1/2)), x1, w1 ), x2, w2 ) ) );
F7 = matlabFunction( simplify( fourier( fourier( (heaviside(x1+1.5) - heaviside(x1+0.5))*(heaviside(x2+1.5) - heaviside(x2+0.5)), x1, w1 ), x2, w2 ) ) );
F8 = matlabFunction( simplify( fourier( fourier( (heaviside(x1+1/5) - heaviside(x1-1/5))*(heaviside(x2+4/5) - heaviside(x2+3/5)), x1, w1 ), x2, w2 ) ) );
xvec6 = F6( multi_idx(:,1) + eps, multi_idx(:,2) + eps );
xvec7 = F7( multi_idx(:,1) + eps, multi_idx(:,2) + eps );
xvec8 = F8( multi_idx(:,1) + eps, multi_idx(:,2) + eps );
% xvec = ones( size( multi_idx( :, 1 ) ) );
% xvec( isnan( xvec ) ) = 0; xvec( isinf( xvec ) ) = 1e10;

% maybe worth looking at the expm_conditioning of the result to ensure 
% that it looks okay, also there are faster and stabler ways to compute 
% this (balance and then exponentiate...
expmXp = expm( Xnum * T ); 
helper6 = expmXp * xvec6;
helper7 = expmXp * xvec7;
helper8 = expmXp * xvec8;
NSamples = 101;
x = linspace( -pi, pi, NSamples );
[ u, v ] = meshgrid( x );
mesaout6 = zeros( NSamples, NSamples );
mesaout7 = zeros( NSamples, NSamples );
mesaout8 = zeros( NSamples, NSamples );
for n = 1:NSamples
    for l = 1:NSamples
%         for k = 1:length( multi_idx )
%             mesaout( n, l ) = mesaout( n, l ) + exp( 1i * u( n, l ) * multi_idx( k, 1 ) ) * exp( 1i * v( n, l ) * multi_idx( k, 2 ) ) * helper( k )/(4*pi^2);
%         end
        mesaout6( n, l ) = (4*pi^2)^-1 * sum( helper6 .* ( exp( 1i * u( n, l ) * multi_idx( :, 1 ) ) .* exp( 1i * v( n, l ) * multi_idx( :, 2 ) ) ) );
        mesaout7( n, l ) = (4*pi^2)^-1 * sum( helper7 .* ( exp( 1i * u( n, l ) * multi_idx( :, 1 ) ) .* exp( 1i * v( n, l ) * multi_idx( :, 2 ) ) ) );
        mesaout8( n, l ) = (4*pi^2)^-1 * sum( helper8 .* ( exp( 1i * u( n, l ) * multi_idx( :, 1 ) ) .* exp( 1i * v( n, l ) * multi_idx( :, 2 ) ) ) );
    
    end
end
figure; hold on;
contour( u, v, real( mesaout6 ) ); colorbar;
xlabel('x1');
ylabel('x2');
shading interp;
plot( Y( :, 1 ), Y( :, 2 ), 'g' );
figure; hold on;
contour( u, v, real( mesaout7 ) ); colorbar;
xlabel('x1');
ylabel('x2');
shading interp;
plot( Y( :, 1 ), Y( :, 2 ), 'g' );
figure; hold on;
contour( u, v, real( mesaout8 ) ); colorbar;
xlabel('x1');
ylabel('x2');
shading interp;
% plot( Y( :, 1 ), Y( :, 2 ), 'g' );
% figure;
% surf( u, v, imag( mesaout ) ); colorbar;
% xlabel('x1');
% ylabel('x2');

%% plot AND SAVE mesa function evolution
% F6 = matlabFunction( simplify( fourier( fourier( (heaviside(x1+3) - heaviside(x1+2))*(heaviside(x2+1/2) - heaviside(x2-1/2)), x1, w1 ), x2, w2 ) ) );
% F7 = matlabFunction( simplify( fourier( fourier( (heaviside(x1+3) - heaviside(x1+2))*(heaviside(x2+3) - heaviside(x2+2)), x1, w1 ), x2, w2 ) ) );
% F8 = matlabFunction( simplify( fourier( fourier( (heaviside(x1+1/2) - heaviside(x1-1/2))*(heaviside(x2+1/2) - heaviside(x2-1/2)), x1, w1 ), x2, w2 ) ) );
% xvec6 = F6( multi_idx(:,1) + eps, multi_idx(:,2) + eps );
% xvec7 = F7( multi_idx(:,1) + eps, multi_idx(:,2) + eps );
% xvec8 = F8( multi_idx(:,1) + eps, multi_idx(:,2) + eps );
% % xvec = ones( size( multi_idx( :, 1 ) ) );
% % xvec( isnan( xvec ) ) = 0; xvec( isinf( xvec ) ) = 1e10;
% 
% % maybe worth looking at the expm_conditioning of the result to ensure 
% % that it looks okay, also there are faster and stabler ways to compute 
% % this (balance and then exponentiate...
% T = linspace( 0, -0.75, 25 );
% for k = 1:length( T )
%     expmXp = expm( Xnum * T( k ) );
%     helper6 = expmXp * xvec6;
%     helper7 = expmXp * xvec7;
%     helper8 = expmXp * xvec8;
%     NSamples = 201;
%     x = linspace( -pi, pi, NSamples );
%     [ u, v ] = meshgrid( x );
%     mesaout6 = zeros( NSamples, NSamples );
%     mesaout7 = zeros( NSamples, NSamples );
%     mesaout8 = zeros( NSamples, NSamples );
%     for n = 1:NSamples
%         for l = 1:NSamples
%             %         for k = 1:length( multi_idx )
%             %             mesaout( n, l ) = mesaout( n, l ) + exp( 1i * u( n, l ) * multi_idx( k, 1 ) ) * exp( 1i * v( n, l ) * multi_idx( k, 2 ) ) * helper( k )/(4*pi^2);
%             %         end
%             mesaout6( n, l ) = (4*pi^2)^-1 * sum( helper6 .* ( exp( 1i * u( n, l ) * multi_idx( :, 1 ) ) .* exp( 1i * v( n, l ) * multi_idx( :, 2 ) ) ) );
%             mesaout7( n, l ) = (4*pi^2)^-1 * sum( helper7 .* ( exp( 1i * u( n, l ) * multi_idx( :, 1 ) ) .* exp( 1i * v( n, l ) * multi_idx( :, 2 ) ) ) );
%             mesaout8( n, l ) = (4*pi^2)^-1 * sum( helper8 .* ( exp( 1i * u( n, l ) * multi_idx( :, 1 ) ) .* exp( 1i * v( n, l ) * multi_idx( :, 2 ) ) ) );
%             
%         end
%     end
%     h = figure( 1 );
%     contour( u, v, real( mesaout6 ) ); colorbar;
%     xlabel('x1');
%     ylabel('x2');
%     shading interp;
%     filename = sprintf( 'indicator_mesa_images/edge/image_%2.3f.pdf', T( k ) );
%     print( h, '-dpdf', filename );
%     h = figure( 2 );
%     contour( u, v, real( mesaout7 ) ); colorbar;
%     xlabel('x1');
%     ylabel('x2');
%     shading interp;
%     filename = sprintf( 'indicator_mesa_images/corner/image_%2.3f.pdf', T( k ) );
%     print( h, '-dpdf', filename );
%     h = figure( 3 );
%     contour( u, v, real( mesaout8 ) ); colorbar;
%     xlabel('x1');
%     ylabel('x2');
%     shading interp;
%     filename = sprintf( 'indicator_mesa_images/center/image_%2.3f.pdf', T( k ) );
%     print( h, '-dpdf', filename );
% end
% % figure;
% % surf( u, v, imag( mesaout ) ); colorbar;
% % xlabel('x1');
% % ylabel('x2');
% 
% 
