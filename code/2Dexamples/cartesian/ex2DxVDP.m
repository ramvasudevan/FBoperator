%% parameters
NBases = 15;
compute_nulls = 0;
scaling = 3;

%% generate multi-idx to make clear what each element in row is
[ gridx, gridy ] = meshgrid( -NBases:NBases );
multi_idx = zeros( size( gridx, 1 ) * size( gridx, 2 ), 2 );
multi_idx( :, 1 ) = gridx( : );
multi_idx( :, 2 ) = gridy( : );

%% compute f_{1-5} and its DFT
% constructing the signals
tic
syms x1 x2 w1 w2;
rot = pi;
h1 = scaling*x2*(heaviside(x1 - pi) - heaviside(x1 + pi))*(heaviside(x2 - pi) - heaviside(x2 + pi));
h2 = -(scaling*x1 + scaling*x2*((scaling*x1)^2 - 1))*(heaviside(x1 - pi) - heaviside(x1 + pi))*(heaviside(x2 - pi) - heaviside(x2 + pi));
% h1 = scaling*x2*exp((-x1^2-x2^2)*5);
% h2 = -(scaling*x1 + scaling*x2*((scaling*x1)^2 - 1))*exp((-x1^2-x2^2)*5);
% h1 = scaling*x2*cos(x1*x2/4);
% h2 = -(scaling*x1 + scaling*x2*((scaling*x1)^2 - 1))*cos(x1*x2/4);
d1 = diff( h1, x1 );
d2 = diff( h2, x2 );
f1 = matlabFunction( scaling * h1, 'vars', { x1; x2 } );
f2 = matlabFunction( scaling * h2, 'vars', { x1; x2 } );
F1 = matlabFunction( simplify( fourier( fourier( -h1, x1, w1 ), x2, w2 ) ), 'vars', { w1; w2} );
F2 = matlabFunction( simplify( fourier( fourier( -d1, x1, w1 ), x2, w2 ) ), 'vars', { w1; w2} );
F3 = matlabFunction( simplify( fourier( fourier( -h2, x1, w1 ), x2, w2 ) ), 'vars', { w1; w2} );
F4 = matlabFunction( simplify( fourier( fourier( -d2, x1, w1 ), x2, w2 ) ), 'vars', { w1; w2} );
f1vec = F1( gridx + eps, gridy + eps );
f2vec = F2( gridx + eps, gridy + eps );
f3vec = F3( gridx + eps, gridy + eps );
f4vec = F4( gridx + eps, gridy + eps );
% f1vec( isnan( f1vec ) ) = 0; f1vec( isinf( f1vec ) ) = 1e10;
% f2vec( isnan( f2vec ) ) = 0; f2vec( isinf( f2vec ) ) = 1e10;
% % f3vec( isnan( f3vec ) ) = 0; f3vec( isinf( f3vec ) ) = 1e10;
% f4vec( isnan( f4vec ) ) = 0; f4vec( isinf( f4vec ) ) = 1e10;
toc
%% compute the eigenrepresentation of X
tic
X = zeros( size( multi_idx, 1 ) );
Z = zeros( size( multi_idx, 1 ) );
s_idx = zeros( size( gridx, 1 ) * size( gridx, 2 ), 2 );
for k = 1:size( multi_idx, 1 )
    c_idx = multi_idx( k, : );
    s_idx( :, 1 ) = c_idx( 1 ) + gridx( : );
    s_idx( :, 2 ) = c_idx( 2 ) + gridy( : );
    f_idx = intersect( intersect( find( c_idx( 1 ) + gridx( : ) >= -NBases ), find( c_idx( 1 ) + gridx( : ) <= NBases ) ), ...
        intersect( find( c_idx( 2 ) + gridy( : ) >= -NBases ), find( c_idx( 2 ) + gridy( : )  <= NBases ) ) );
    h_idx = any( all( bsxfun( @eq, reshape( s_idx.', 1, 2, [] ), multi_idx), 2 ), 3 );
    X( h_idx, k ) = 1i * c_idx( 1 ) * f1vec( f_idx ) + f2vec( f_idx ) + ...
        1i * c_idx( 2 ) * f3vec( f_idx ) + f4vec( f_idx );
end
toc
%% compute eigenvalues of different portions
if ( compute_nulls )
    [ helperL, helperR ] = spspaces( X, 3 );
%     [ ~, helper ] = spspaces( Z, 2 );
    VR = helperR{1};
    zeroeigR = helperR{3};
    VL = helperL{1};
    zeroeigL = helperL{3};
else
    [ VR, D ] = eig( X, eye( size( X ) ) );
    
    zeroeigR = find( abs( diag( D ) ) <= 1e-10 );
end

%% just for plotting
tic
% NSamples = 3;
% x = linspace( -1, 1, NSamples );
% [ u, v ] = meshgrid( x );
% helper = reshape( exp( 1i * kron( u( : ), multi_idx( :, 1 ) ) ).* exp( 1i * kron( v( : ), multi_idx( :, 2 ) ) ), ...
%     size( multi_idx, 1 ), NSamples, NSamples );
% out = zeros( length( zeroeig ), NSamples, NSamples );
% for m = 1:length( zeroeig )
%     for n = 1:NSamples
%         for l = 1:NSamples0
%             out( m, n, l ) = V( :, zeroeig( m ) )' * helper( :, l, n );
%         end
%     end
% end
NSamples = 101;
x = linspace( -pi, pi, NSamples );
[ u, v ] = meshgrid( x );
z_idx = any( all( bsxfun( @eq, reshape( [ 0 0 ].', 1, 2, [] ), multi_idx), 2 ), 3 );
outR = zeros( length( zeroeigR ), NSamples, NSamples );
% outL = zeros( length( zeroeigL ), NSamples, NSamples );
for m = 1:length( zeroeigR )
    for n = 1:NSamples
        for l = 1:NSamples
%             out( m, n, l ) = ( exp( 1i * u( n, l ) * multi_idx( :, 1 ) ) .* exp( 1i * v( n, l ) * multi_idx( :, 2 ) ) )' * V( :, zeroeig( m ) )/V( z_idx, zeroeig( m ) );
            outR( m, n, l ) = (4*pi^2)^-1 * sum( VR( :, zeroeigR( m ) ) .* ( exp( 1i * u( n, l ) * multi_idx( :, 1 ) ) .* exp( 1i * v( n, l ) * multi_idx( :, 2 ) ) ) ) ;
%             outL( m, n, l ) = VL( zeroeigL( m ), : ) * ( exp( 1i * u( n, l ) * multi_idx( :, 1 ) ) .* exp( 1i * v( n, l ) * multi_idx( :, 2 ) ) );
%             out( m, n, l ) = ( exp( 1i * u( n, l ) * multi_idx( :, 1 ) ) .* exp( 1i * v( n, l ) * multi_idx( :, 2 ) ) )' * V( zeroeig( m ), : )';
        
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
    pause;
end

%% plot vector field
[ u, v ] = meshgrid( -pi:0.05*pi:( pi ), -pi:0.05*pi:( pi ) );
Du = f1( u, v );
Dv = f2( u, v );
figure; hold on;
quiver( u, v, Du, Dv, 5 )
[ ~, Y ] = ode45( @vdpvf, [0 1000], [0.75 0.5] );
plot( Y( :, 1 ), Y( :, 2 ), 'g' );

%% plot mesa function evolution
% T = 1e-5;
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
% expmXp = expm( X * T ); 
% helper6 = expmXp * xvec6;
% helper7 = expmXp * xvec7;
% helper8 = expmXp * xvec8;
% NSamples = 201;
% x = linspace( -pi, pi, NSamples );
% [ u, v ] = meshgrid( x );
% mesaout6 = zeros( NSamples, NSamples );
% mesaout7 = zeros( NSamples, NSamples );
% mesaout8 = zeros( NSamples, NSamples );
% for n = 1:NSamples
%     for l = 1:NSamples
% %         for k = 1:length( multi_idx )
% %             mesaout( n, l ) = mesaout( n, l ) + exp( 1i * u( n, l ) * multi_idx( k, 1 ) ) * exp( 1i * v( n, l ) * multi_idx( k, 2 ) ) * helper( k )/(4*pi^2);
% %         end
%         mesaout6( n, l ) = (4*pi^2)^-1 * sum( helper6 .* ( exp( 1i * u( n, l ) * multi_idx( :, 1 ) ) .* exp( 1i * v( n, l ) * multi_idx( :, 2 ) ) ) );
%         mesaout7( n, l ) = (4*pi^2)^-1 * sum( helper7 .* ( exp( 1i * u( n, l ) * multi_idx( :, 1 ) ) .* exp( 1i * v( n, l ) * multi_idx( :, 2 ) ) ) );
%         mesaout8( n, l ) = (4*pi^2)^-1 * sum( helper8 .* ( exp( 1i * u( n, l ) * multi_idx( :, 1 ) ) .* exp( 1i * v( n, l ) * multi_idx( :, 2 ) ) ) );
%     
%     end
% end
% figure;
% contour( u, v, real( mesaout6 ) ); colorbar;
% xlabel('x1');
% ylabel('x2');
% shading interp;
% figure;
% contour( u, v, real( mesaout7 ) ); colorbar;
% xlabel('x1');
% ylabel('x2');
% shading interp;
% figure;
% contour( u, v, real( mesaout8 ) ); colorbar;
% xlabel('x1');
% ylabel('x2');
% shading interp;
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
% T = linspace( 0, 0.75, 25 );
% for k = 1:length( T )
%     expmXp = expm( X * T( k ) );
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
