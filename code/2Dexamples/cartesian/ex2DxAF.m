%% parameters
NBases = 11;
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
% h1 = -x1 * (heaviside( x1 + pi ) - heaviside( x1 - pi )) * (heaviside( x2 + pi ) - heaviside( x2 - pi ));
% h2 = -x2 * (heaviside( x1 + pi ) - heaviside( x1 - pi )) * (heaviside( x2 + pi ) - heaviside( x2 - pi ));
% h1 = -x1 * taylor( exp(-1/(1-(x1/pi)^100)) ) * taylor( exp(-1/(1-(x2/pi)^100)) ) * (exp(1))^2;
% h2 = -x2 * taylor( exp(-1/(1-(x1/pi)^100)) ) * taylor( exp(-1/(1-(x2/pi)^100)) ) * (exp(1))^2;
d1 = diff( -x1, x1 );
d2 = diff( -x2, x2 );
% h3 = heaviside( pi^2/4 - ( x1^2 + x2^2 ) );
% f1 = matlabFunction( -scaling * x1 * h1 * h2 );
% f2 = matlabFunction( -scaling * x2 * h1 * h2 );
% F1 = matlabFunction( simplify( fourier( fourier( scaling * x1 * h1 * h2, x1, w1 ), x2, w2 ) ) );
% F2 = matlabFunction( fourier( fourier( scaling * h1 * h2, x1, w1 ), x2, w2 ) );
% F3 = matlabFunction( fourier( fourier( scaling * x1 * h2 * d1, x1, w1 ), x2, w2 ) );
% F4 = matlabFunction( simplify( fourier( fourier( scaling * x2 * h1 * h2, x1, w1 ), x2, w2 ) ) );
% F5 = matlabFunction( fourier( fourier( scaling * h1 * h2, x1, w1 ), x2, w2 ) );
% F6 = matlabFunction( fourier( fourier( scaling * x2 * h1 * d2, x1, w1 ), x2, w2 ) );
% f1vec = F1( gridx + eps, gridy + eps );
% f2vec = F2( gridx + eps, gridy + eps );
% f3vec = F3( gridx + eps, gridy + eps );
% f4vec = F4( gridx + eps, gridy + eps );
% f5vec = F5( gridx + eps, gridy + eps );
% f6vec = F6( gridx + eps, gridy + eps );
f1 = matlabFunction( scaling * h1 );
f2 = matlabFunction( scaling * h2 );
F1 = matlabFunction( simplify( fourier( fourier( scaling * h1, x1, w1 ), x2, w2 ) ), 'vars', {w1,w2} );
F2 = matlabFunction( simplify( fourier( fourier( scaling * d1, x1, w1 ), x2, w2 ) ), 'vars', {w1,w2} );
F3 = matlabFunction( simplify( fourier( fourier( scaling * h2, x1, w1 ), x2, w2 ) ), 'vars', {w1,w2} );
F4 = matlabFunction( simplify( fourier( fourier( scaling * d2, x1, w1 ), x2, w2 ) ), 'vars', {w1,w2} );
% F1 = matlabFunction( simplify( int( int( scaling * h1 * h2 * exp(-1i*w1*x1)*exp(-1i*w2*x2), x1, -pi,pi ), x2, -pi,pi ) ) );
% F2 = matlabFunction( simplify( int( int( scaling * d1 * h2 * exp(-1i*w1*x1)*exp(-1i*w2*x2), x1, -pi,pi ), x2, -pi,pi ) ) );
% F3 = matlabFunction( simplify( int( int( scaling * h1 * h2 * exp(-1i*w1*x1)*exp(-1i*w2*x2), x1, -pi,pi ), x2, -pi,pi ) ) );
% F4 = matlabFunction( simplify( int( int( scaling * h1 * d2 * exp(-1i*w1*x1)*exp(-1i*w2*x2), x1, -pi,pi ), x2, -pi,pi ) ) );
f1vec = F1( gridx + eps, gridy + eps );
f2vec = F2( gridx + eps, gridy + eps );
f3vec = F3( gridx + eps, gridy + eps );
f4vec = F4( gridx + eps, gridy + eps );
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
%     X( h_idx, k ) = 1i * c_idx( 1 ) * f1vec( f_idx ) + f2vec( f_idx ) + f3vec( f_idx ) + ...
%         1i * c_idx( 2 ) * f4vec( f_idx ) + f5vec( f_idx ) + f6vec( f_idx );
    X( h_idx, k ) = 1i * c_idx( 1 ) * f1vec( f_idx ) + f2vec( f_idx ) + ...
        1i * c_idx( 2 ) * f3vec( f_idx ) + f4vec( f_idx );
%     for j = 1:size( multi_idx, 1 )
%         r = multi_idx( j, 1 );
%         s = multi_idx( j, 2 );
%         y = multi_idx( j, 1 ) - c_idx( 1 );
%         z = multi_idx( j, 2 ) - c_idx( 2 );
%         if( y >= -NBases && y <= NBases && z >= -NBases && z <= NBases )
%             Z( j, k ) = 4 * 1i * sin( y + eps ) * sin( z + eps )/( ( y + eps ) * ( z + eps ) ) * ( c_idx( 1 ) + c_idx( 2 ) + y + z );
%         end
%     end
end
toc
%% compute eigenvalues of different portions
if ( compute_nulls )
    [ ~, helper ] = spspaces( X, 2 );
%     [ ~, helper ] = spspaces( Z, 2 );
    V = helper{1};
    zeroeig = helper{3};
else
    [ V, D ] = eig( X, eye( size( X ) ) );
    
    zeroeig = find( abs( diag( D ) ) <= 1e-10 );
    
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
out = zeros( length( zeroeig ), NSamples, NSamples );
for m = 1:length( zeroeig )
    for n = 1:NSamples
        for l = 1:NSamples
%             out( m, n, l ) = ( exp( 1i * u( n, l ) * multi_idx( :, 1 ) ) .* exp( 1i * v( n, l ) * multi_idx( :, 2 ) ) )' * V( :, zeroeig( m ) )/V( z_idx, zeroeig( m ) );
            out( m, n, l ) = ( exp( 1i * u( n, l ) * multi_idx( :, 1 ) ) .* exp( 1i * v( n, l ) * multi_idx( :, 2 ) ) )' * V( :, zeroeig( m ) );
            
            %             out( m, n, l ) = ( exp( 1i * u( n, l ) * multi_idx( :, 1 ) ) .* exp( 1i * v( n, l ) * multi_idx( :, 2 ) ) )' * f1vec( : );
            
        end
    end
end
toc
%% actual plotting
for m = 1:length( zeroeig )
    % for m = 75
    figure;
    surf( u, v, real( squeeze( out( m, :, : ) ) ) ); colorbar;
    xlabel('x1');
    ylabel('x2');
    shading interp;
    figure;
    contour( u, v, imag( squeeze( out( m, :, : ) ) ) ); colorbar;
    xlabel('x1');
    ylabel('x2');
    shading interp;
    %     pause;
end

%% plot vector field
[ u, v ] = meshgrid( -pi:0.1*pi:( pi ), -pi:0.1*pi:( pi ) );
Du = zeros( size( u ) );
Dv = zeros( size( v ) );
% for j = 1:size( u, 1 )
%     for k = 1:size( u, 2 )
%         Du( j, k ) = f3( u( j, k ), v( j, k ) );
%         Dv( j, k ) = f3( u( j, k ), v( j, k ) );
% 
%     end
% end
Du = f1( u, v );
Dv = f2( u, v );
figure; quiver( u, v, Du, Dv )

