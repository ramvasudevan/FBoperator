%% parameters
NBases = 12;
compute_nulls = 1;
scaling = 20;

%% generate multi-idx to make clear what each element in row is
[ gridx, gridy ] = meshgrid( -NBases:NBases );
multi_idx = zeros( size( gridx, 1 ) * size( gridx, 2 ), 2 );
multi_idx( :, 1 ) = gridx( : );
multi_idx( :, 2 ) = gridy( : );

%% compute f_{1-5} and its DFT
% constructing the signals
tic
syms x1 x2 w1 w2;
h1 = heaviside( x1 + pi/2 ) - heaviside( x1 - pi/2 );
h2 = heaviside( x2 + pi/2 ) - heaviside( x2 - pi/2 );
d1 = dirac( x1 + pi/2 ) - dirac( x1 - pi/2 );
d2 = dirac( x2 + pi/2 ) - dirac( x2 - pi/2 );
f1 = matlabFunction( fourier( fourier( scaling * x2 * h1 * h2, x1, w1 ), x2, w2 ) );
f2 = matlabFunction( fourier( fourier( scaling * x2 * d1 * h2, x1, w1 ), x2, w2 ) );
f3 = matlabFunction( fourier( fourier( scaling * x1 * h1 * h2 + (scaling * x1)^2 * scaling * x2 * h1 * h2 - ...
    scaling * x2 * h1 * h2, x1, w1 ), x2, w2 ) );
f4 = matlabFunction( fourier( fourier( scaling * x1 * h1 * d2 + (scaling * x1)^2 * scaling * x2 * h1 * d2 - ...
    scaling * x2 * h1 * d2, x1, w1 ), x2, w2 ) );
f5 = matlabFunction( fourier( fourier( (scaling * x1)^2 * scaling * h1 * h2 - scaling * h1 * h2, x1, w1 ), x2, w2 ) );
f1vec = f1( gridx + eps, gridy + eps );
f2vec = f2( gridx + eps, gridy + eps );
f3vec = f3( gridx + eps, gridy + eps );
f4vec = f4( gridx + eps, gridy + eps );
f5vec = f5( gridx + eps, gridy + eps );
toc
%% compute the eigenrepresentation of X
tic
X = zeros( size( multi_idx, 1 ) );
s_idx = zeros( size( gridx, 1 ) * size( gridx, 2 ), 2 );
for k = 1:size( multi_idx, 1 )
    c_idx = multi_idx( k, : );
    s_idx( :, 1 ) = c_idx( 1 ) + gridx( : );
    s_idx( :, 2 ) = c_idx( 2 ) + gridy( : );
    f_idx = intersect( intersect( find( c_idx( 1 ) + gridx( : ) >= -NBases ), find( c_idx( 1 ) + gridx( : ) <= NBases ) ), ...
        intersect( find( c_idx( 2 ) + gridy( : ) >= -NBases ), find( c_idx( 2 ) + gridy( : )  <= NBases ) ) );
    h_idx = any( all( bsxfun( @eq, reshape( s_idx.', 1, 2, [] ), multi_idx), 2 ), 3 );
    X( h_idx, k ) = 1i * c_idx( 1 ) * f1vec( f_idx ) + f2vec( f_idx ) - 1i * c_idx( 2 ) * f3vec( f_idx ) ...
        - f4vec( f_idx ) - f5vec( f_idx );
end
toc
%% compute eigenvalues of different portions
if ( compute_nulls )
    [ ~, helper ] = spspaces( X, 2 );
    V = helper{1};
    zeroeig = helper{3};
else
    [ V, D ] = eig( X, eye( size( X ) ) );
    
    zeroeig = [];
    for k = 1:length( D )
        if ( D( k, k ) == 0 )
            zeroeig = [ zeroeig; k ];
        end
    end
end

%% just for plotting
tic
NSamples = 101;
x = linspace( -pi, pi, NSamples );
[ u, v ] = meshgrid( x );
out = zeros( length( zeroeig ), NSamples, NSamples );
for m = 1:length( zeroeig )
    for n = 1:NSamples
        for l = 1:NSamples
            out( m, n, l ) = ( exp( 1i * u( n, l ) * multi_idx( :, 1 ) ) .* exp( 1i * v( n, l ) * multi_idx( :, 2 ) ) )' * V( :, zeroeig( m ) );
            %             out( m, n, l ) = ( exp( 1i * u( n, l ) * multi_idx( :, 1 ) ) .* exp( 1i * v( n, l ) * multi_idx( :, 2 ) ) )' * f1vec( : );
            
        end
    end
end
toc
%% actual plotting
for m = 1:size( zeroeig, 2 )
    % for m = 75
    figure;
    mesh( x, x, real( squeeze( out( m, :, : ) ) ) ); colorbar;
    %     pause;
end


