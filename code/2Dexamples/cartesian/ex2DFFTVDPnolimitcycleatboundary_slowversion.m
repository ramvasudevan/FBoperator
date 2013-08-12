%% parameters
NBases = 31;
compute_nulls = 0;
NSamples = 101;
scaling = 1.25;
bumpx = 2.9;
bumpy = 2.9;
expbump = 200;
vdpscale = 0.5;

%% generate multi-idx to make clear what each element in row is
[ gridx, gridy ] = meshgrid( -NBases:NBases );
multi_idx = zeros( size( gridx, 1 ) * size( gridx, 2 ), 2 );
multi_idx( :, 1 ) = gridx( : );
multi_idx( :, 2 ) = gridy( : );

%% compute f_{1-5} and its DFT
% constructing the signals
syms x1 x2 e p
tic
Fs = 4096; % increase this number to increase the condition number of the matrix (i.e. increase the accuracy of our computation)
t = -pi:( ( 2 * pi )/Fs ):( pi - ( 2 * pi)/Fs );
% t = 0:( ( 2 * pi )/Fs ):( 2 * pi - ( 2 * pi)/Fs );
[ tx, ty ] = meshgrid( t );
% ORIGINAL VDP
% a = scaling * x2;
% b = -(scaling * x1 + scaling * x2 * ( ( scaling * x1 )^2 - 1 ));
% % h1 = matlabFunction( ( 1/exp(1/(1-(x1/p)^100)) * ( a + ( e^-1 - 1/exp(1/(1-(x2/p)^100) ) ) * b ) * e ), 'vars', {x1,x2,e,p} );
% % h2 = matlabFunction( ( 1./exp(1/(1-(x2./p).^100)) .* ( b + ( e^-1 - 1./exp(1./(1-(x1./p).^100) ) ) .* a ) * e ), 'vars', {x1,x2,e,p}  );
% % d1 = matlabFunction( ( diff( 1/exp(1./(1-(x1./p).^100)) .* ( a + ( e^-1 - 1./exp(1./(1-(x2/p).^100) ) ) .* b ) * e, x1 ) ), 'vars', {x1,x2,e,p} );
% % d2 = matlabFunction( ( diff( 1/exp(1./(1-(x2./p).^100)) .* ( b + ( e^-1 - 1./exp(1./(1-(x1/p).^100) ) ) .* a ) * e, x2 ) ), 'vars', {x1,x2,e,p} );
% h1 = @(x1,x2) exp(1).*exp(-1./(1-(x1./pi).^100)).*(3.*x2 + (3.*x1 + 3.*x2.*(9.*x1.^2 - 1)).*(exp(-1./(1-(x2./pi).^100)) - exp(-1)));
% h2 = @(x1,x2) -exp(1).*exp(-1./(1-(x2./pi).^100)).*(3.*x1 + 3.*x2.*(9.*x1.^2 - 1) + 3.*x2.*(exp(-1./(1-(x1./pi).^100)) - exp(-1)));
% d1 = @(x1,x2) exp(1).*exp(-1./(1-(x1./pi).^100)).*(exp(-1./(1-(x2./pi).^100)) - exp(-1)).*(54.*x1.*x2 + 3) - (100.*exp(1).*x1.^99.*exp(-1./(1-(x1./pi).^100)).*(3.*x2 + (3.*x1 + 3.*x2.*(9.*x1.^2 - 1)).*(exp(-1./(1-(x2./pi).^100)) - exp(-1))))./(pi^100.*(x1.^100./pi^100 - 1).^2);
% d2 = @(x1,x2) (100.*exp(1).*x2.^99.*exp(-1./(1-(x2./pi).^100)).*(3.*x1 + 3.*x2.*(9.*x1.^2 - 1) + 3.*x2.*(exp(-1./(1-(x1./pi).^100)) - exp(-1))))./(pi^100.*(x2.^100./pi^100 - 1).^2) - exp(1).*exp(-1./(1-(x2./pi).^100)).*(3.*exp(-1./(1-(x1./pi).^100)) - 3*exp(-1) + 27.*x1.^2 - 3);
% DOUBLE GYRE
% a = pi * sin( x1 ) .* cos( x2 );
% b = -pi * cos( x1 ) .* sin( x2 );
% h1 = @(x1,x2) exp(1)*exp(-1./(1-(x1./pi).^100)).*(pi.*cos(x2).*sin(x1) + pi.*cos(x1).*sin(x2).*(exp(-1./(1-(x2./pi).^100)) - exp(-1)));
% h2 = @(x1,x2) -exp(1)*exp(-1./(1-(x2./pi).^100)).*(pi.*cos(x1).*sin(x2) + pi.*cos(x2).*sin(x1).*(exp(-1./(1-(x1./pi).^100)) - exp(-1)));
% d1 = @(x1,x2) exp(1)*exp(-1./(1-(x1./pi).^100)).*(pi.*cos(x1).*cos(x2) - pi.*sin(x1).*sin(x2).*(exp(-1./(1-(x2./pi).^100)) - exp(-1))) - (100.*exp(1).*x1.^99.*exp(-1./(1-(x1./pi).^100)).*(pi.*cos(x2).*sin(x1) + pi.*cos(x1).*sin(x2).*(exp(-1./(1-(x2./pi).^100)) - exp(-1))))./(pi^100.*(x1.^100/pi^100 - 1).^2);
% d2 = @(x1,x2) (100.*exp(1).*x2.^99.*exp(-1./(1-(x2./pi).^100)).*(pi.*cos(x1).*sin(x2) + pi.*cos(x2).*sin(x1).*(exp(-1./(1-(x1./pi).^100)) - exp(-1))))./(pi^100.*(x2.^100/pi^100 - 1).^2) - exp(1).*exp(-1./(1-(x2./pi).^100)).*(pi.*cos(x1).*cos(x2) - pi.*sin(x1).*sin(x2).*(exp(-1./(1-(x1./pi).^100)) - exp(-1)));
% SCALED VDP (limit cycle around 2)
% a = scaling * x2;
% b = -(scaling * x1 + 0.1 * scaling * x2 * ( ( scaling * x1 )^2 - 1 ));
% h1 = @(x1,x2) exp(1)*exp(-1./(1-(x1./pi).^100)).*(3.*x2 + (3.*x1 + (3.*x2.*(9.*x1.^2 - 1))/10).*(exp(-1./(1-(x2./pi).^100)) - exp(-1)));
% h2 = @(x1,x2) -exp(1)*exp(-1./(1-(x2./pi).^100)).*(3.*x1 + (3.*x2.*(9.*x1.^2 - 1))/10 + 3.*x2.*(exp(-1./(1-(x1./pi).^100)) - exp(-1)));
% d1 = @(x1,x2) exp(1)*exp(-1./(1-(x1./pi).^100)).*(exp(-1./(1-(x2./pi).^100)) - exp(-1)).*((27.*x1.*x2)/5 + 3) - (100.*exp(1).*x1.^99.*exp(-1./(1-(x1./pi).^100)).*(3.*x2 + (3.*x1 + (3.*x2.*(9.*x1.^2 - 1))/10).*(exp(-1./(1-(x2./pi).^100)) - exp(-1))))./(pi^100*(x1.^100/pi^100 - 1).^2);
% d2 = @(x1,x2) (100*exp(1).*x2.^99.*exp(-1./(1-(x2./pi).^100)).*(3.*x1 + (3.*x2.*(9.*x1.^2 - 1))/10 + 3.*x2.*(exp(-1./(1-(x1./pi).^100)) - exp(-1))))./(pi^100*(x2.^100/pi^100 - 1).^2) - exp(1).*exp(-1./(1-(x2./pi).^100)).*(3*exp(-1./(1-(x1./pi).^100)) - 3*exp(-1) + (27*x1.^2)/10 - 3/10);
% SCALED VDP smooth bump function in each direction
% a = scaling * x2;
% b = -(scaling * x1 + 0.1 * scaling * x2 * ( ( scaling * x1 )^2 - 1 ));
% % h1 = matlabFunction( ( 1/exp(1/(1-(x1/p)^100)) * 1/exp(1/(1-(x2/p)^100)) * a * e^2 ), 'vars', {x1,x2,e,p} );
% % h2 = matlabFunction( ( 1/exp(1/(1-(x1/p)^100)) * 1/exp(1/(1-(x2/p)^100)) * b * e^2 ), 'vars', {x1,x2,e,p}  );
% % d1 = matlabFunction( ( diff( 1/exp(1/(1-(x1/p)^100)) * 1/exp(1/(1-(x2/p)^100)) * a * e^2, x1 ) ), 'vars', {x1,x2,e,p} );
% % d2 = matlabFunction( ( diff( 1/exp(1/(1-(x1/p)^100)) * 1/exp(1/(1-(x2/p)^100)) * b * e^2, x2 ) ), 'vars', {x1,x2,e,p} );
% h1 = @(x1,x2) exp(1)^2.*x2.*exp(-1./(1-(x1./pi).^100)).*exp(-1./(1-(x2./pi).^100));
% h2 = @(x1,x2) -exp(1)^2.*exp(-1./(1-(x1./pi).^100)).*exp(-1./(1-(x2./pi).^100)).*(x1 + (x2.*(x1.^2 - 1))/10);
% d1 = @(x1,x2) -(100*exp(1)^2.*x1.^99.*x2.*exp(-1./(1-(x1./pi).^100)).*exp(-1./(1-(x2./pi).^100)))./(pi^100.*(x1.^100/pi^100 - 1).^2);
% d2 = @(x1,x2) (100*exp(1)^2.*x2.^99.*exp(-1./(1-(x1./pi).^100)).*exp(-1./(1-(x2./pi).^100)).*(x1 + (x2.*(x1.^2 - 1))/10))./(pi^100.*(x2.^100/pi^100 - 1).^2) - exp(1)^2.*exp(-1./(1-(x1./pi).^100)).*exp(-1./(1-(x2./pi).^100)).*(x1.^2/10 - 1/10);
% Original VDP smooth bump function in each direction
% a = scaling * x2;
% b = -(scaling * x1 + scaling * x2 * ( ( scaling * x1 )^2 - 1 ));
% h1 = @(x1,x2) 2.*exp(1).^2.*x2.*exp(-1./(1-(x1./pi).^100)).*exp(-1./(1-(x2./pi).^100));
% h2 = @(x1,x2) -exp(1).^2.*exp(-1./(1-(x1./pi).^100)).*exp(-1./(1-(x2./pi).^100)).*(2.*x1 + 2.*x2.*(4.*x1.^2 - 1));
% d1 = @(x1,x2) -(200.*exp(1).^2.*x1.^99.*x2.*exp(-1./(1-(x1./pi).^100)).*exp(-1./(1-(x2./pi).^100)))./(pi.^100.*(x1.^100./pi.^100 - 1).^2);
% d2 = @(x1,x2) (100.*exp(1).^2.*x2.^99.*exp(-1./(1-(x1./pi).^100)).*exp(-1./(1-(x2./pi).^100)).*(2.*x1 + 2.*x2.*(4.*x1.^2 - 1)))./(pi.^100.*(x2.^100./pi.^100 - 1).^2) - exp(1).^2.*exp(-1./(1-(x1./pi).^100)).*exp(-1./(1-(x2./pi).^100)).*(8.*x1.^2 - 2);
% original VDP different bump function in each direction
% h1 = matlabFunction( scaling * x2 * exp( -(x1/bumpx)^expbump ) * exp( -(x2/bumpy)^expbump ) );
% h2 = matlabFunction( -(scaling*x1 + scaling*vdpscale*x2.*((scaling*x1).^2 - 1)) * exp( -(x1/bumpx)^expbump ) * exp( -(x2/bumpy)^expbump ) );
% d1 = matlabFunction( diff( scaling * x2 * exp( -(x1/bumpx)^expbump ) * exp( -(x2/bumpy)^expbump ), x1 ) );
% d2 = matlabFunction( diff( -(scaling*x1 + scaling*vdpscale*x2.*((scaling*x1).^2 - 1)) * exp( -(x1/bumpx)^expbump ) * exp( -(x2/bumpy)^expbump ), x2 ) );
h1 = matlabFunction( scaling * x2 * ( pi^2 - x2^2 ) * ( pi^2 - x1^2 ) );
h2 = matlabFunction( -(scaling*x1 + scaling*vdpscale*x2.*((scaling*x1).^2 - 1)) * ( pi^2 - x2^2 ) * ( pi^2 - x1^2 ) );
d1 = matlabFunction( diff( scaling * x2 * ( pi^2 - x2^2 ) * ( pi^2 - x1^2 ), x1 ) );
d2 = matlabFunction( diff( -(scaling*x1 + scaling*vdpscale*x2.*((scaling*x1).^2 - 1)) * ( pi^2 - x2^2 ) * ( pi^2 - x1^2 ), x2 ) );
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
F1 = fftshift( fft2( fftshift( f1 ), nfft, nfft ) )/nfft^2 * (2 * pi);
F2 = fftshift( fft2( fftshift( f2 ), nfft, nfft ) )/nfft^2 * (2 * pi);
F3 = fftshift( fft2( fftshift( f3 ), nfft, nfft ) )/nfft^2 * (2 * pi);
F4 = fftshift( fft2( fftshift( f4 ), nfft, nfft ) )/nfft^2 * (2 * pi);
% F1 = ( fft2( f1, nfft, nfft ) )/nfft^2;
% F2 = ( fft2( f2, nfft, nfft ) )/nfft^2;
% F3 = ( fft2( f3, nfft, nfft ) )/nfft^2;
% F4 = ( fft2( f4, nfft, nfft ) )/nfft^2;
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
for k = 1:size( multi_idx, 1 )
    s_idx = zeros( size( gridx, 1 ) * size( gridx, 2 ), 2 );
    c_idx = multi_idx( k, : );
    s_idx( :, 1 ) = c_idx( 1 ) + gridx( : );
    s_idx( :, 2 ) = c_idx( 2 ) + gridy( : );
    f_idx = intersect( intersect( find( c_idx( 1 ) + gridx( : ) >= -NBases ), find( c_idx( 1 ) + gridx( : ) <= NBases ) ), ...
        intersect( find( c_idx( 2 ) + gridy( : ) >= -NBases ), find( c_idx( 2 ) + gridy( : )  <= NBases ) ) );
    h_idx = any( all( bsxfun( @eq, reshape( s_idx.', 1, 2, [] ), multi_idx), 2 ), 3 );
    %     X( h_idx, k ) = 1i * c_idx( 1 ) * f1vec( f_idx ) + f2vec( f_idx ) + f3vec( f_idx ) + ...
    %         1i * c_idx( 2 ) * f4vec( f_idx ) + f5vec( f_idx ) + f6vec( f_idx );
    % correct non-parallellizable version
    %     Xnum( h_idx, k ) = 1i * c_idx( 1 ) * f1vec( f_idx ) + f2vec( f_idx ) + ...
    %         1i * c_idx( 2 ) * f3vec( f_idx ) + f4vec( f_idx );
    blah = zeros( size( multi_idx, 1 ), 1 );
    blah( h_idx ) = 1i * c_idx( 1 ) * f1vec( f_idx ) + f2vec( f_idx ) + ...
        1i * c_idx( 2 ) * f3vec( f_idx ) + f4vec( f_idx );
    Xnum( :, k ) = blah;
end
toc
%% compute eigenvalues of different portions
if ( compute_nulls )
    [ helperL, helperR ] = spspaces( Xnum, 3, 1 );
    VR = helperR{1};
    zeroeigR = helperR{3}
    VL = helperL{1};
    zeroeigL = helperL{3};
else
    tic
    [ VR, D ] = eig( Xnum );
    
    zeroeigR = find( abs( real( diag( D ) ) ) < 1e-3 )
    toc
%     tic;
%     [UR, D, VR ] = svd( Xnum );
%     zeroeigR = find( abs( real( diag( D ) ) ) < 1e-1 )
%     toc
end

%% just for plotting
tic
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

%% plot vector field
% [ u, v ] = meshgrid( -pi:0.05*pi:( pi ), -pi:0.05*pi:( pi ) );
% Du = zeros( size( u ) );
% Dv = zeros( size( v ) );
% for j = 1:size( u, 1 )
%     for k = 1:size( u, 2 )
%         blah = vdpvf( 0, [ u( j, k ); v( j, k ) ] );
%         Du( j, k ) = blah( 1 );
%         Dv( j, k ) = blah( 2 );
%     end
% end
Du = h1( u, v );
Dv = h2( u, v );
figure; hold on;
quiver( u, v, Du, Dv );
[ ~, Y ] = ode45( @vdpvf, [0 10], [ -0.1 0 ] );
plot( Y( :, 1 ), Y( :, 2 ), 'g' );
% print( h, '-dpdf', 'indicator_mesa_images/vf_and_lc.pdf' );
% [ ~, Y ] = ode45( @vdpvf, [0 1000], [ 0  pi-1.5e-1 ] );
% plot( Y( :, 1 ), Y( :, 2 ), 'r' );

%% actual plotting
x = linspace( -pi, pi, NSamples );
[ u, v ] = meshgrid( x );
for m = 1:length( zeroeigR )
    % for m = 75
    figure; hold on;
    surf( u, v, real( squeeze( outR( m, :, : ) ) ) ); colorbar;
    plot( Y( :, 1 ), Y( :, 2 ), 'k', 'LineWidth', 2 );
    xlabel('x1');
    ylabel('x2');
    shading interp;
    %     figure; hold on;
    %     surf( u, v, imag( squeeze( outR( m, :, : ) ) ) ); colorbar;
    %         plot( Y( :, 1 ), Y( :, 2 ), 'k', 'LineWidth', 2 );
    %     xlabel('x1');
    %     ylabel('x2');
    %     shading interp;
end

%% plot mesa function evolution
syms x1 x2 w1 w2;
annout = 2.1;
annin = 1.9;
T = -1;
% F6 = matlabFunction( simplify( fourier( fourier( (heaviside(x1+1.5) - heaviside(x1+0.5))*(heaviside(x2+1/2) - heaviside(x2-1/2)), x1, w1 ), x2, w2 ) ) );
F6 = matlabFunction( simplify( fourier( fourier( (heaviside(x1+3) - heaviside(x1+2))*(heaviside(x2+1/2) - heaviside(x2-1/2)), x1, w1 ), x2, w2 ) ) );
F7 = matlabFunction( simplify( fourier( fourier( (heaviside(x1+1.5) - heaviside(x1+0.5))*(heaviside(x2+1.5) - heaviside(x2+0.5)), x1, w1 ), x2, w2 ) ) );
% F8 = matlabFunction( simplify( fourier( fourier( (heaviside(x1+4/5) - heaviside(x1+2/5) + heaviside(x1-2/5) - heaviside(x1-4/5))*(heaviside(x2+4/5) - heaviside(x2+2/5) + heaviside(x2-2/5) - heaviside(x2-4/5)), x1, w1 ), x2, w2 ) ) );
xvec6 = F6( multi_idx(:,1) + eps, multi_idx(:,2) + eps );
xvec7 = F7( multi_idx(:,1) + eps, multi_idx(:,2) + eps );
% xvec8 = F8( multi_idx(:,1) + eps, multi_idx(:,2) + eps );

h8 = @(x1,x2) 50*(heaviside(annout^2-x1.^2-x2.^2) - heaviside(annin^2-x1.^2-x2.^2));
f8 = h8( tx, ty );
F8 = fftshift( fft2( fftshift( f8 ), nfft, nfft ) )/nfft^2;
xvec8 = interp2( freqsx, freqsy, F8, gridx, gridy, 'nearest' );
xvec8 = reshape( xvec8, size( xvec8, 1 ) * size( xvec8, 2 ), 1 );
% xvec = ones( size( multi_idx( :, 1 ) ) );
% xvec( isnan( xvec ) ) = 0; xvec( isinf( xvec ) ) = 1e10;

% maybe worth looking at the expm_conditioning of the result to ensure
% that it looks okay, also there are faster and stabler ways to compute
% this (balance and then exponentiate...
expmXp = expm( Xnum * T );
helper6 = expmXp * xvec6;
helper7 = expmXp * xvec7;
helper8 = expmXp * xvec8;
% NSamples = 201;
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
% figure; hold on;
% contour( u, v, real( mesaout6 ) ); colorbar;
% xlabel('x1');
% ylabel('x2');
% shading interp;
% plot( Y( :, 1 ), Y( :, 2 ), 'k', 'LineWidth', 2 );
% figure; hold on;
% contour( u, v, real( mesaout7 ) ); colorbar;
% xlabel('x1');
% ylabel('x2');
% shading interp;
% plot( Y( :, 1 ), Y( :, 2 ), 'k', 'LineWidth', 2 );
figure; hold on;
surf( u, v, real( mesaout8 ) ); colorbar;
xlabel('x1');
ylabel('x2');
shading interp;
plot( Y( :, 1 ), Y( :, 2 ), 'k', 'LineWidth', 2 );

%% determine if BRS is actually returning the correct thing
th = 0.6;
getthere = zeros( size( u ) );
for j = 1:size( u, 1 )
    for k = 1:size( u, 2 )
        [ ~, P ] = ode45( @vdpvf, [T 0], [ u( j, k ) v( j, k ) ] );
        getthere( j, k ) = sqrt( P( end, 1 )^2 + P( end, 2 )^2 );
    end
end
length( find( mesaout8 > th ) )
length( intersect( find( getthere > annin ), intersect( find( getthere < annout ), find( mesaout8 > th ) ) ) )


% th = 0.25;
% [ indx, indy ] = find( mesaout8 > th );
% getthere = zeros( size( indx ) );
% for k = 1:length( indx )
%     [ ~, P ] = ode45( @vdpvf, [T 0], [ u( indx( k ), indy( k ) ) v( indx( k ), indy( k ) ) ] );
%     %     getthere( k ) = h8( P( end, 1 ), P( end, 2 ) );
%     getthere( k ) = sqrt( P( end, 1 )^2 + P( end, 2 )^2 );
% end

%% plot AND SAVE mesa function evolution
% syms x1 x2 w1 w2;
% F6 = matlabFunction( simplify( fourier( fourier( (heaviside(x1+3) - heaviside(x1+2))*(heaviside(x2+1/2) - heaviside(x2-1/2)), x1, w1 ), x2, w2 ) ) );
% F7 = matlabFunction( simplify( fourier( fourier( (heaviside(x1+3) - heaviside(x1+2))*(heaviside(x2+3) - heaviside(x2+2)), x1, w1 ), x2, w2 ) ) );
% % F8 = matlabFunction( simplify( fourier( fourier( (heaviside(x1+1/2) - heaviside(x1-1/2))*(heaviside(x2+1/2) - heaviside(x2-1/2)), x1, w1 ), x2, w2 ) ) );
% xvec6 = F6( multi_idx(:,1) + eps, multi_idx(:,2) + eps );
% xvec7 = F7( multi_idx(:,1) + eps, multi_idx(:,2) + eps );
% % xvec8 = F8( multi_idx(:,1) + eps, multi_idx(:,2) + eps );
%
% h8 = @(x1,x2) 10*(heaviside(1.5^2-x1.^2-x2.^2) - heaviside(1^2-x1.^2-x2.^2));
% f8 = h8( tx, ty );
% F8 = fftshift( fft2( fftshift( f8 ), nfft, nfft ) )/nfft^2;
% xvec8 = interp2( freqsx, freqsy, F8, gridx, gridy, 'nearest' );
% xvec8 = reshape( xvec8, size( xvec8, 1 ) * size( xvec8, 2 ), 1 );
% % xvec = ones( size( multi_idx( :, 1 ) ) );
% % xvec( isnan( xvec ) ) = 0; xvec( isinf( xvec ) ) = 1e10;
%
% % maybe worth looking at the expm_conditioning of the result to ensure
% % that it looks okay, also there are faster and stabler ways to compute
% % this (balance and then exponentiate...
% T = linspace( 0, 5, 100 );
% for k = 1:length( T )
%     expmXp = expm( Xnum * T( k ) );
% %     helper6 = expmXp * xvec6;
% %     helper7 = expmXp * xvec7;
%     helper8 = expmXp * xvec8;
%     NSamples = 201;
%     x = linspace( -pi, pi, NSamples );
%     [ u, v ] = meshgrid( x );
% %     mesaout6 = zeros( NSamples, NSamples );
% %     mesaout7 = zeros( NSamples, NSamples );
%     mesaout8 = zeros( NSamples, NSamples );
%     for n = 1:NSamples
%         for l = 1:NSamples
%             %         for k = 1:length( multi_idx )%% parameters
