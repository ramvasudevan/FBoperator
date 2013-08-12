%% parameters
NBases = 6;
compute_nulls = 1;
scaling = 1;

%% generate multi-idx to make clear what each element in row is
[ gridx, gridy ] = meshgrid( -NBases:NBases );
multi_idx = zeros( size( gridx, 1 ) * size( gridx, 2 ), 2 );
multi_idx( :, 1 ) = gridx( : );
multi_idx( :, 2 ) = gridy( : );

%% compute f_{1-5} and its DFT
% constructing the signals
tic
syms x1 x2 w1 w2;
% h1 = ( heaviside( x1 + 1 ) - heaviside( x1 - 1 ) ) * ( heaviside( x2 + 1 ) - heaviside( x2 - 1 ) );
% h2 = ( heaviside( x1 + 1 ) - heaviside( x1 - 1 ) ) * ( heaviside( x2 + 1 ) - heaviside( x2 - 1 ) );
% h1 = -x1 * exp((-x1^6-x2^6)*5);
% h2 = -x2 * exp((-x1^6-x2^6)*5);
% h1 = -x1 * ( heaviside( x1 + 1 - x2 ) - heaviside( x1 - 1 - x2 ) ) * ( heaviside( -x1 + 1 - x2 ) - heaviside( -x1 - 1 - x2 ) );
% h2 = -x2 * ( heaviside( x1 + 1 - x2 ) - heaviside( x1 - 1 - x2 ) ) * ( heaviside( -x1 + 1 - x2 ) - heaviside( -x1 - 1 - x2 ) );
% h1 = -x1 * exp((-x1^2-x2^2)*5) * cos( pi/4 ) + x2 * exp((-x1^2-x2^2)*5) * sin( pi/4 );
% h2 = -x1 * exp((-x1^2-x2^2)*5) * sin( pi/4 ) - x2 * exp((-x1^2-x2^2)*5) * cos( pi/4 );
h1 = -x1 * exp(-1/(1-(x1/pi)^2)) * exp(-1/(1-(x2/pi)^2)) * (exp(1))^2;
h2 = -x2 * exp(-1/(1-(x1/pi)^2)) * exp(-1/(1-(x2/pi)^2)) * (exp(1))^2;
d1 = diff( h1 );
d2 = diff( h2 );
% f1 = matlabFunction( scaling * h1, 'vars', { x1; x2 } );
% f2 = matlabFunction( scaling * h2, 'vars', { x1; x2 } );
f1 = @(x1,x2) -x1 * exp(-1/(1-(x1/pi)^2)) * exp(-1/(1-(x2/pi)^2)) * (exp(1))^2;
f2 = @(x1,x2) -x2 * exp(-1/(1-(x1/pi)^2)) * exp(-1/(1-(x2/pi)^2)) * (exp(1))^2;
% f3 = matlabFunction( h3, 'vars', { x1, x2 } );
% h11 = ( x1 * exp( -x1^2/ 5 ) )^2;
% h12 = ( x1 * exp( -x2^2/ 5 ) )^2;
% h21 = ( x2 * exp( -x1^2/ 5 ) )^2;
% h22 = ( x2 * exp( -x2^2/ 5 ) )^2;
% % h1 = ( heaviside( 1 - x1^2 - x2^2 ) );
% % h2 = ( heaviside( 1 - x1^2 - x2^2 ) );
% d11 = 2*x1*exp(-(2*x1^2)/5) - (4*x1^3*exp(-(2*x1^2)/5))/5;
% d12 = 2*x1*exp(-(2*x2^2)/5);
% d21 = 2*x2*exp(-(2*x1^2)/5);
% d22 = 2*x2*exp(-(2*x2^2)/5) - (4*x2^3*exp(-(2*x2^2)/5))/5;
% f1 = @(x1,x2) -min( ( x1 * exp( -x1^2/ 5 ) )^2, ( x1 * exp( -x2^2/ 5 ) )^2 );
% f2 = @(x1,x2) -min( ( x2 * exp( -x1^2/ 5 ) )^2, ( x2 * exp( -x2^2/ 5 ) )^2 );

% numerically computing the fourier coefficients instead of symbolically
% computing them! NOTE DOES NOT WORK IF DELTA FUNCTIONS START APPEARING in
% derivatives!!
[ numx, numy ] = meshgrid( linspace( -pi, pi, 1001 ) );
dxdy = ( numx( 1, 2 ) - numx( 1, 1 ) )^2;
g1vec = zeros( size( gridx ) );
g2vec = zeros( size( gridx ) );
g3vec = zeros( size( gridx ) );
g4vec = zeros( size( gridx ) );
for j = 1:size( gridx, 1 )
    tic
    for k = 1:size( gridx, 2 )
%         F1 = matlabFunction( simplify( -scaling * h1 * exp( -1i * gridx( j, k ) * x1 ) * exp( -1i * gridy( j, k ) * x2 ) ) );
%         F2 = matlabFunction( simplify( -scaling * d1 * exp( -1i * gridx( j, k ) * x1 ) * exp( -1i * gridy( j, k ) * x2 ) ) );
%         F3 = matlabFunction( simplify( -scaling * h2 * exp( -1i * gridx( j, k ) * x1 ) * exp( -1i * gridy( j, k ) * x2 ) ) );
%         F4 = matlabFunction( simplify( -scaling * d2 * exp( -1i * gridx( j, k ) * x1 ) * exp( -1i * gridy( j, k ) * x2 ) ) );
%         g1vec( j, k ) = sum( sum( F1( numx, numy ) ) ) * dxdy;
%         g2vec( j, k ) = sum( sum( F2( numx, numy ) ) ) * dxdy;
%         g3vec( j, k ) = sum( sum( F3( numx, numy ) ) ) * dxdy;
%         g4vec( j, k ) = sum( sum( F4( numx, numy ) ) ) * dxdy;
        F1 = @(x1,x2)  -x1 .* exp(-1./(1-(x1/pi)^2)) .* exp(-1./(1-(x2/pi)^2)) .* (exp(1))^2.* exp( -1i * gridx( j, k ) * x1 ) .* exp( -1i * gridy( j, k ) * x2 );
        F2 = @(x1,x2) (exp(-1./(1-x1^2/pi^2)).*(x1^4-4*pi^2*x1^2+pi^4))./(pi^2-x1^2)^2.* exp(-1./(1-(x2/pi)^2)).* exp( -1i * gridx( j, k ) * x1 ) .* exp( -1i * gridy( j, k ) * x2 );
        F3 = @(x1,x2) -x2 .* exp(-1./(1-(x1/pi)^2)) .* exp(-1./(1-(x2/pi)^2)) .* (exp(1))^2.* exp( -1i * gridx( j, k ) * x1 ) .* exp( -1i * gridy( j, k ) * x2 );
        F4 = @(x1,x2) (exp(-1./(1-x2^2/pi^2)) .* (x2^4-4*pi^2*x2^2+pi^4))./(pi^2-x2^2)^2 .* exp(-1./(1-(x1/pi)^2)).* exp( -1i * gridx( j, k ) * x1 ) .* exp( -1i * gridy( j, k ) * x2 );
        g1vec( j, k ) = sum( sum( F1( numx, numy ) ) ) * dxdy;
        g2vec( j, k ) = sum( sum( F2( numx, numy ) ) ) * dxdy;
        g3vec( j, k ) = sum( sum( F3( numx, numy ) ) ) * dxdy;
        g4vec( j, k ) = sum( sum( F4( numx, numy ) ) ) * dxdy;
    end
    toc
    j
end
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
    Xnum( h_idx, k ) = 1i * c_idx( 1 ) * g1vec( f_idx ) + g2vec( f_idx ) + ...
        1i * c_idx( 2 ) * g3vec( f_idx ) + g4vec( f_idx );
%     for j = 1:size( multi_idx, 1 )
%         r = multi_idx( j, 1 );
%         s = multi_idx( j, 2 );
%         y = multi_idx( j, 1 ) - c_idx( 1 );
%         z = multi_idx( j, 2 ) - c_idx( 2 );
%         if( y >= -NBases && y <= NBases && z >= -NBases && z <= NBases )
%             Z( j, k ) = 1i * c_idx( 1 ) * F1( y + eps, z + eps ) + F2( y + eps, z + eps ) ...
%                 + 1i * c_idx( 2 ) * F3( y + eps, z + eps ) + F4( y + eps, z + eps );
%         end
%     end
end
toc
%% compute eigenvalues of different portions
if ( compute_nulls )
    [ helperL, helperR ] = spspaces( Xnum, 3, 1 );
%     [ ~, helper ] = spspaces( Z, 2 );
    VR = helperR{1};
    zeroeigR = helperR{3};
    VL = helperL{1};
    zeroeigL = helperL{3};
else
    [ V, D ] = eig( X, eye( size( Xnum ) ) );
    
    zeroeig = [];
    for k = 1:length( D )
        if ( D( k, k ) == 0 )
            zeroeig = [ zeroeig; k ];
        end
    end
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
%         for l = 1:NSamples
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
%     figure;
%     surf( u, v, real( squeeze( outL( m, :, : ) ) ) ); colorbar;
%     xlabel('x1');
%     ylabel('x2');
%     shading interp;
%     figure;
%     contour( u, v, imag( squeeze( outL( m, :, : ) ) ) ); colorbar;
%     xlabel('x1');
%     ylabel('x2');
%     shading interp;
    %     pause;
end

%% plot vector field
[ u, v ] = meshgrid( -pi:0.05*pi:( pi ), -pi:0.05*pi:( pi ) );
Du = zeros( size( u ) );
Dv = zeros( size( v ) );
for j = 1:size( u, 1 )
    for k = 1:size( u, 2 )
        Du( j, k ) = f1( u( j, k ), v( j, k ) );
        Dv( j, k ) = f2( u( j, k ), v( j, k ) );

    end
end
figure; quiver( u, v, Du, Dv )

%% plot mesa function evolution
% T = 1e-3;
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
% expmXp = expm( Xnum * T ); 
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
% surf( u, v, real( mesaout6 ) ); colorbar;
% xlabel('x1');
% ylabel('x2');
% shading interp;
% figure;
% surf( u, v, real( mesaout7 ) ); colorbar;
% xlabel('x1');
% ylabel('x2');
% shading interp;
% figure;
% surf( u, v, real( mesaout8 ) ); colorbar;
% xlabel('x1');
% ylabel('x2');
% shading interp;
% % figure;
% % surf( u, v, imag( mesaout ) ); colorbar;
% % xlabel('x1');
% % ylabel('x2');

