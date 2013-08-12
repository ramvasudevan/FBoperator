%% parameters
NBases = 10;
compute_nulls = 0;
NSamples = 101;
% scaling = 1;
% vdpscale = 1;
addpath('ex_vf');

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
[ tx, ty ] = meshgrid( t );
[ f1, f3, f2, f4 ] = vdpR2( tx, ty );
% computing their FFTs
nfft = length( t );
F1 = fftshift( fft2( fftshift( -f1 ), nfft, nfft ) )/nfft^2 * (4*pi^2);
F2 = fftshift( fft2( fftshift( -f2 ), nfft, nfft ) )/nfft^2 * (4*pi^2);
F3 = fftshift( fft2( fftshift( -f3 ), nfft, nfft ) )/nfft^2 * (4*pi^2);
F4 = fftshift( fft2( fftshift( -f4 ), nfft, nfft ) )/nfft^2 * (4*pi^2);
% frequency vector
freqs = ( -nfft/2:( nfft/2 - 1 ) ) * Fs/nfft ;
sampled_freqs = -NBases:NBases;
[ freqsx, freqsy ] = meshgrid( freqs );
fvec = zeros( 2 * NBases + 1, 2 * NBases + 1, 4 );
fvec( :, :, 1 ) = interp2( freqsx, freqsy, F1, gridx, gridy, 'nearest' );
fvec( :, :, 2 ) = interp2( freqsx, freqsy, F2, gridx, gridy, 'nearest' );
fvec( :, :, 3 ) = interp2( freqsx, freqsy, F3, gridx, gridy, 'nearest' );
fvec( :, :, 4 ) = interp2( freqsx, freqsy, F4, gridx, gridy, 'nearest' );
gvec = zeros( NBases + 1, 2 * NBases + 1, 4 );
for k = 1:4
        gvec( 1:NBases, :, k ) = 1i * ( fvec(  ( 2 * NBases + 2 ) - (1:NBases), :, k ) - fvec( 1:NBases, :, k ) );
%     gvec( 1:NBases, :, k ) = ( fvec(  ( 2 * NBases + 2 ) - (1:NBases), :, k ) + fvec( 1:NBases, :, k ) );
    gvec( NBases + 1, :, k ) = 2 * fvec( NBases + 1, :, k );
end
hvec = zeros( NBases + 1, NBases + 1, 4 );
for k = 1:4
    hvec( :, 1:NBases, k ) = 1i * ( gvec(  :, ( 2 * NBases + 2 ) - (1:NBases), k ) - gvec( :, 1:NBases, k ) );
%     hvec( :, 1:NBases, k ) = ( gvec(  :, ( 2 * NBases + 2 ) - (1:NBases), k ) + gvec( :, 1:NBases, k ) );
    hvec( :, NBases + 1, k ) = 2 * gvec( :, NBases + 1, k );
end
toc
%% compute the new basis that we are working with
B = zeros( ( 2 * NBases + 1 )^2, ( 2 * NBases + 1 )^2 );
helper = zeros( ( 2 * NBases + 1 )^2, 1 );
helper( [ 1 2 ( 2 * NBases + 1 ) + 1 ( 2 * NBases + 1 ) + 2 ] ) = 1;
for k = 1:length( multi_idx )
    B( k:length( multi_idx ), k ) = helper( 1:( length(k:length( multi_idx ) ) ) );
end
A = inv( B );
%% compute the eigenrepresentation of X
tic
Xnum = zeros( size( multi_idx, 1 ) );
for k = 1:size( multi_idx, 1 )
    s_idx = zeros( size( gridx, 1 ) * size( gridx, 2 ), 2 );
    c_idx = multi_idx( k, : );
    s_idx( :, 1 ) = c_idx( 1 ) + gridx( : );
    s_idx( :, 2 ) = c_idx( 2 ) + gridy( : );
    f_idx = s_idx( :, 1 ) >= -NBases & s_idx( :, 1 ) <= NBases & s_idx( :, 2 ) >= -NBases & s_idx( :, 2 ) <= NBases;
    helper = ( s_idx( :, 1 ) + NBases ) * ( 2 * NBases + 1 ) + ( s_idx( :, 2 ) + NBases + 1 );
    h_idx = helper( helper > 0 & f_idx );
    
    blah = zeros( size( multi_idx, 1 ), 1 );
    blah( h_idx ) = 1i * c_idx( 1 ) * f1vec( f_idx ) + f2vec( f_idx ) + ...
        1i * c_idx( 2 ) * f3vec( f_idx ) + f4vec( f_idx );
    Xnum( :, k ) = A * blah;
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
    %     tic
    %     [ VR, D ] = eig( Xnum );
    %     zeroeigR = find( abs( real( diag( D ) ) ) < 1e-3 )
    %     toc
    %     tic;
    %     [UR, D, VR ] = svd( Xnum );
    %     zeroeigR = find( abs( real( diag( D ) ) ) < 1e-1 )
    %     toc

    tic
    zeroeigR = 1:5;
    [ VR, D ] = eigs( Xnum, length( zeroeigR ), 'SM' );
    eigens = diag( D );
    toc;
end

%% just for plotting
% tic
% x = linspace( -pi, pi , NSamples );
% [ u, v ] = meshgrid( x );
% outR = zeros( length( zeroeigR ), NSamples, NSamples );
% % outL = zeros( length( zeroeigL ), NSamples, NSamples );
% parfor m = 1:length( zeroeigR )
%     for n = 1:NSamples
%         for l = 1:NSamples
%             for m = 1:length( h )
%             outC( n, l ) = outC( n, l ) + h( m ) * sum( exp( 1i * u( n, l ) * multi_idx( find( B( :, m  ) ), 1  ) ) .* exp( 1i * v( n, l ) * multi_idx( find( B( :, m  ) ), 2  ) ) );
%         end
%         outC( n, l ) = outC( n, l ) * 0.25 * (4*pi^2)^-1;
%             outR( m, n, l ) = (4*pi^2)^-1 * sum( VR( :, zeroeigR( m ) ) .* ( exp( 1i * u( n, l ) * multi_idx( :, 1 ) ) .* exp( 1i * v( n, l ) * multi_idx( :, 2 ) ) ) ) ;
%         end
%     end
% end
% toc

%% plot vector field
[ Du, Dv, ~, ~ ] = vdpR2( u, v );
figure; hold on;
quiver( u, v, Du, Dv );
[ ~, Y ] = ode23( @vdpR2_plot, [0 1e-18], [ 0 -0.1 ] );
plot( Y( :, 1 ), Y( :, 2 ), 'g' );
% print( h, '-dpdf', 'indicator_mesa_images/vf_and_lc.pdf' );
% [ ~, Y ] = ode45( @vdpvf, [0 1000], [ 0  pi-1.5e-1 ] );
% plot( Y( :, 1 ), Y( :, 2 ), 'r' );

%% actual plotting
% for m = 1:length( zeroeigR )
%     % for m = 75
%     h = figure; hold on;
%     surf( u, v, real( squeeze( outR( m, :, : ) ) ) ); colorbar;
% %     plot( Y( :, 1 ), Y( :, 2 ), 'k', 'LineWidth', 2 );
%     xlabel('x1');
%     ylabel('x2');
%     shading interp;
% %     filename = sprintf( 'indicator_mesa_images/VDP_30/eigdist_image_%3.10f_%3.10f.pdf', real( eigens( m ) ), imag( eigens( m ) ) );
% %     print( h, '-dpdf', filename );
%     figure; hold on;
%     surf( u, v, imag( squeeze( outR( m, :, : ) ) ) ); colorbar;
%     plot( Y( :, 1 ), Y( :, 2 ), 'k', 'LineWidth', 2 );
%     xlabel('x1');
%     ylabel('x2');
%     shading interp;
% end

%% plot mesa function evolution
% syms x1 x2 w1 w2;
% T = 1e-6;
% F6 = matlabFunction( simplify( fourier( fourier( (heaviside(x1+2.5) - heaviside(x1+1.5))*(heaviside(x2+1/2) - heaviside(x2-1/2)), x1, w1 ), x2, w2 ) ) );
% F7 = matlabFunction( simplify( fourier( fourier( (heaviside(x1+0.5) - heaviside(x1-0.5))*(heaviside(x2+2.5) - heaviside(x2+1.5)), x1, w1 ), x2, w2 ) ) );
% F8 = matlabFunction( simplify( fourier( fourier( (heaviside(x1+1/2) - heaviside(x1-1/2))*(heaviside(x2+1/2) - heaviside(x2-1/2)), x1, w1 ), x2, w2 ) ) );
% xvec6 = F6( multi_idx(:,1) + eps, multi_idx(:,2) + eps );
% xvec7 = F7( multi_idx(:,1) + eps, multi_idx(:,2) + eps );
% xvec8 = F8( multi_idx(:,1) + eps, multi_idx(:,2) + eps );
%
% % maybe worth looking at the expm_conditioning of the result to ensure
% % that it looks okay, also there are faster and stabler ways to compute
% % this (balance and then exponentiate...
% expmXp = expm( Xnum * T );
% helper6 = expmXp * xvec6;
% helper7 = expmXp * xvec7;
% helper8 = expmXp * xvec8;
% mesaout6 = zeros( NSamples, NSamples );
% mesaout7 = zeros( NSamples, NSamples );
% mesaout8 = zeros( NSamples, NSamples );
% for n = 1:NSamples
%     for l = 1:NSamples
%         %         for k = 1:length( multi_idx )
%         %             mesaout( n, l ) = mesaout( n, l ) + exp( 1i * u( n, l ) * multi_idx( k, 1 ) ) * exp( 1i * v( n, l ) * multi_idx( k, 2 ) ) * helper( k )/(4*pi^2);
%         %         end
%         mesaout6( n, l ) = (4*pi^2)^-1 * sum( helper6 .* ( exp( 1i * u( n, l ) * multi_idx( :, 1 ) ) .* exp( 1i * v( n, l ) * multi_idx( :, 2 ) ) ) );
%         mesaout7( n, l ) = (4*pi^2)^-1 * sum( helper7 .* ( exp( 1i * u( n, l ) * multi_idx( :, 1 ) ) .* exp( 1i * v( n, l ) * multi_idx( :, 2 ) ) ) );
%         mesaout8( n, l ) = (4*pi^2)^-1 * sum( helper8 .* ( exp( 1i * u( n, l ) * multi_idx( :, 1 ) ) .* exp( 1i * v( n, l ) * multi_idx( :, 2 ) ) ) );
%
%     end
% end
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
% figure; hold on;
% contour( u, v, real( mesaout8 ) ); colorbar;
% xlabel('x1');
% ylabel('x2');
% shading interp;
% plot( Y( :, 1 ), Y( :, 2 ), 'k', 'LineWidth', 2 );

%% determine if BRS is actually returning the correct thing
% th = 0.6;
% getthere = zeros( size( u ) );
% for j = 1:size( u, 1 )
%     for k = 1:size( u, 2 )
%         [ ~, P ] = ode45( @vdpvf, [T 0], [ u( j, k ) v( j, k ) ] );
%         getthere( j, k ) = sqrt( P( end, 1 )^2 + P( end, 2 )^2 );
%     end
% end
% length( find( mesaout8 > th ) )
% length( intersect( find( getthere > annin ), intersect( find( getthere < annout ), find( mesaout8 > th ) ) ) )


% th = 0.25;
% [ indx, indy ] = find( mesaout8 > th );
% getthere = zeros( size( indx ) );
% for k = 1:length( indx )
%     [ ~, P ] = ode45( @vdpvf, [T 0], [ u( indx( k ), indy( k ) ) v( indx( k ), indy( k ) ) ] );
%     %     getthere( k ) = h8( P( end, 1 ), P( end, 2 ) );
%     getthere( k ) = sqrt( P( end, 1 )^2 + P( end, 2 )^2 );
% end

%% plot AND SAVE mesa function evolution
syms x1 x2 w1 w2;
F6 = matlabFunction( simplify( fourier( fourier( (heaviside(x1+1.75) - heaviside(x1+0.75))*(heaviside(x2+1/2) - heaviside(x2-1/2)), x1, w1 ), x2, w2 ) ) );
% F7 = matlabFunction( simplify( fourier( fourier( (heaviside(x1+0.5) - heaviside(x1-0.5))*(heaviside(x2+2.5) - heaviside(x2+1.5)), x1, w1 ), x2, w2 ) ) );
% F8 = matlabFunction( simplify( fourier( fourier( (heaviside(x1+1/2) - heaviside(x1-1/2))*(heaviside(x2+1/2) - heaviside(x2-1/2)), x1, w1 ), x2, w2 ) ) );
xvec6 = F6( multi_idx(:,1) + eps, multi_idx(:,2) + eps );
% xvec7 = F7( multi_idx(:,1) + eps, multi_idx(:,2) + eps );
% xvec8 = F8( multi_idx(:,1) + eps, multi_idx(:,2) + eps );
% xvec = ones( size( multi_idx( :, 1 ) ) );
% xvec( isnan( xvec ) ) = 0; xvec( isinf( xvec ) ) = 1e10;

% maybe worth looking at the expm_conditioning of the result to ensure
% that it looks okay, also there are faster and stabler ways to compute
% this (balance and then exponentiate...
T = linspace( 0, 1e-19, 20 );
% tic
% mXp = expm( Xnum );
% toc

x = linspace( -pi, pi, NSamples );
[ u, v ] = meshgrid( x );
for k = 1:length( T )
    tic
    expmXp = expm( Xnum  * T( k ) );
    toc
    helper6 = expmXp * xvec6;
    %     helper7 = expmXp * xvec7;
    %     helper8 = expmXp * xvec8;
    mesaout6 = zeros( NSamples, NSamples );
    %     mesaout7 = zeros( NSamples, NSamples );
    %     mesaout8 = zeros( NSamples, NSamples );
    for n = 1:NSamples
        for l = 1:NSamples
            %         for k = 1:length( multi_idx )
            %             mesaout( n, l ) = mesaout( n, l ) + exp( 1i * u( n, l ) * multi_idx( k, 1 ) ) * exp( 1i * v( n, l ) * multi_idx( k, 2 ) ) * helper( k )/(4*pi^2);
            %         end
            mesaout6( n, l ) = (4*pi^2)^-1 * sum( helper6 .* ( exp( 1i * u( n, l ) * multi_idx( :, 1 ) ) .* exp( 1i * v( n, l ) * multi_idx( :, 2 ) ) ) );
            %             mesaout7( n, l ) = (4*pi^2)^-1 * sum( helper7 .* ( exp( 1i * u( n, l ) * multi_idx( :, 1 ) ) .* exp( 1i * v( n, l ) * multi_idx( :, 2 ) ) ) );
            %             mesaout8( n, l ) = (4*pi^2)^-1 * sum( helper8 .* ( exp( 1i * u( n, l ) * multi_idx( :, 1 ) ) .* exp( 1i * v( n, l ) * multi_idx( :, 2 ) ) ) );
            
        end
    end
    h = figure( k ); hold on;
    surf( u, v, real( mesaout6 ) ); colorbar;
    xlabel('x1');
    ylabel('x2');
    shading interp;
    plot( Y( :, 1 ), Y( :, 2 ), 'k', 'LineWidth', 2 );
    T( k )
    %     filename = sprintf( 'indicator_mesa_images/scaledVDP_30/image_%2.6f.pdf', T( k ) )
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
end
% figure;
% surf( u, v, imag( mesaout ) ); colorbar;
% xlabel('x1');
% ylabel('x2');
