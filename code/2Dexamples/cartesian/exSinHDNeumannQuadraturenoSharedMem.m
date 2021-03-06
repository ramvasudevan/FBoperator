%% parameters
NBases = 2;
compute_nulls = 0;
NSamplesInt = 100;
NSamplesPlot = 201;
addpath('ex_vf');

%% generate multi-idx to make clear what each element in row is
[ gridx, gridy ] = meshgrid( 0:NBases );
multi_idx = zeros( size( gridx, 1 ) * size( gridx, 2 ), 2 );
multi_idx( :, 1 ) = gridx( : );
multi_idx( :, 2 ) = gridy( : );

%% generate integration parameters
[ t, Weights ] = lgwt( NSamplesInt, -pi, pi );
Weights = Weights * Weights';
Weights = Weights( : );
[ tx, ty ] = meshgrid( t );
tx = tx( : );
ty = ty( : );
%% compute functions for later use
tic
[ f1, f2, f3, f4 ] = SinR2( tx, ty );
Trig = zeros( size( multi_idx, 1 ), NSamplesInt^2 );
DfTrig = zeros( size( multi_idx, 1 ), NSamplesInt^2 );
for k = 1:size( multi_idx, 1 )
    cx = multi_idx( k, 1 );
    cy = multi_idx( k, 2 );
    if( mod( cx, 2 ) && mod( cy, 2 ) )
        Trig( k, : ) = sin( cx * tx/2 ) .* sin( cy * ty/2 );
        DfTrig( k, : ) =  Weights .* ( f1 .* cx/2 .* cos( cx * tx/2 ) .* sin( cy * ty/2 ) + f2 .* cy/2 .* sin( cx * tx/2 ) .* cos( cy * ty/2 ) + ...
            1/2 .* ( f3 + f4 ) .* sin( cx * tx/2 ) .* sin( cy * ty/2 ) );
    elseif( mod( cx, 2 ) )
        Trig( k, : ) = sin( cx * tx/2 ) .* cos( cy * ty/2 );
        DfTrig( k, : ) =  Weights .* ( f1 .* cx/2 .* cos( cx * tx/2 ) .* cos( cy * ty/2 ) + f2 .* -cy/2 .* sin( cx * tx/2 ) .* sin( cy * ty/2 ) + ...
            1/2 .* ( f3 + f4 ) .* sin( cx * tx/2 ) .* cos( cy * ty/2 ) );
    elseif( mod( cy, 2 ) )
        Trig( k, : ) = cos( cx * tx/2 ) .* sin( cy * ty/2 );
        DfTrig( k, : ) = Weights .* ( f1 .* -cx/2 .* sin( cx * tx/2 ) .* sin( cy * ty/2 ) + f2 .* cy/2 .* cos( cx * tx/2 ) .* cos( cy * ty/2 ) + ...
            1/2 .* ( f3 + f4 ) .* cos( cx * tx/2 ) .* sin( cy * ty/2 ) );
    else
        Trig( k, : ) = cos( cx * tx/2 ) .* cos( cy * ty/2 );
        DfTrig( k, : ) = Weights .* ( f1 .* -cx/2 .* sin( cx * tx/2 ) .* cos( cy * ty/2 ) + f2 .* -cy/2 .* cos( cx * tx/2 ) .* sin( cy * ty/2 ) + ...
            1/2 .* ( f3 + f4 ) .* cos( cx * tx/2 ) .* cos( cy * ty/2 ) );
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
X = zeros( size( multi_idx, 1 ) );
for k = 1:size( multi_idx, 1 )
    X( k, : ) = sum( bsxfun( @times, squeeze( DfTrig( k, : ) ), Trig ), 2 )';
end
% Y = X;
% for j = 1:size( multi_idx, 1 )
%     for k = 1:size( multi_idx, 1 )
%         if( multi_idx( j, 1 ) <= NBases - 1 || multi_idx( j, 2 ) <= NBases - 1 || multi_idx( k, 1 ) <= NBases - 1 || multi_idx( k, 2 ) <= NBases - 1 )
%            Y( j, k ) = 0;
%         end
%     end
% end
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
    zeroeig = find( abs( ( eigens ) ) <= 1e-14 )
    toc
    %             tic
    %             [ U, D, V ] = svd( X );
    %             zeroeig = find( abs( real( diag( D ) ) ) < 1e-10 )
    %             toc
    %     tic
    %     zeroeig = 1:10;
    %     [ V, D ] = eigs( X, length( zeroeig ), 'SM' );
    %     eigens = diag( D )
    %     toc;
end

%% just for plotting
tic
x = linspace( -pi, pi, NSamplesPlot );
scaling_factor = ones( size( multi_idx, 1 ), 1 );
% scaling_factor( 1 ) = 1/sqrt( 2 );
[ u, v ] = meshgrid( x );
out = zeros( length( zeroeig ), NSamplesPlot, NSamplesPlot );
for m = 2:2:length( zeroeig )
    for n = 1:NSamplesPlot
        for l = 1:NSamplesPlot
            for k = 1:size( multi_idx, 1 )
                cx = multi_idx( k, 1 );
                cy = multi_idx( k, 2 );
                if( mod( cx, 2 ) && mod( cy, 2 ) )
                    out( m, n, l ) = out( m, n, l ) + V( k, zeroeig( m ) ) * sin( cx * u( n, l )/2 ) .* sin( cy * v( n, l )/2 );
                elseif( mod( cx, 2 ) )
                    out( m, n, l ) = out( m, n, l ) + V( k, zeroeig( m ) ) * sin( cx * u( n, l )/2 ) .* cos( cy * v( n, l )/2 ) * scaling_factor( cy + 1 );
                elseif( mod( cy, 2 ) )
                    out( m, n, l ) = out( m, n, l ) + V( k, zeroeig( m ) ) * cos( cx * u( n, l )/2 ) .* sin( cy * v( n, l )/2 ) * scaling_factor( cx + 1 );
                else
                    out( m, n, l ) = out( m, n, l ) + V( k, zeroeig( m ) ) * cos( cx * u( n, l )/2 ) .* cos( cy * v( n, l )/2 ) * scaling_factor( cx + 1 ) * scaling_factor( cy + 1 );
                end
            end
        end
    end
end
toc

%% plot vector field
[ Du, Dv, ~, ~ ] = SinR2( u, v );
figure; hold on;
quiver( u, v, Du, Dv );
[ ~, Y ] = ode23( @SinR2_plot, [0 0.1], [ -pi/4 -pi/2 ] );
plot( Y( :, 1 ), Y( :, 2 ), 'g' );
% print( h, '-dpdf', 'indicator_mesa_images/vf_and_lc.pdf' );
% [ ~, Y ] = ode45( @vdpvf, [0 1000], [ 0  pi-1.5e-1 ] );
% plot( Y( :, 1 ), Y( :, 2 ), 'r' );

%% actual plotting
for m = 2:2:length( zeroeig )
    % for m = 75
    h = figure; hold on;
    %     set(gcf,'renderer','painters')
    surf( u, v, squeeze( out( m, :, : ) ).* conj( squeeze( out( m, :, : ) ) ) ); colorbar;
    shading flat
    %     contour( u, v, squeeze( out( m, :, : ) ).* conj( squeeze( out( m, :, : ) ) ), 100 ); colorbar;
    plot( Y( :, 1 ), Y( :, 2 ), 'g', 'LineWidth', 10 );
    view( -15, 60 );
    xlabel('x1');
    ylabel('x2');
%     filename = sprintf( 'DG_densities/eigdist_image_%3.25f.pdf', imag( eigens( zeroeig( m ) ) ) );
%     print( h, '-dpdf', filename );
end

%% plot mesa function evolution
syms x1 x2 w1 w2;
T = linspace( 0, 0.1, 10 );
F6 = matlabFunction( simplify( fourier( fourier( (heaviside(x1+1.5) - heaviside(x1+0.5))*(heaviside(x2+1.5) - heaviside(x2+1/2)), x1, w1 ), x2, w2 ) ) );
xvec6 = F6( multi_idx(:,1) + eps, multi_idx(:,2) + eps );

% maybe worth looking at the expm_conditioning of the result to ensure
% that it looks okay, also there are faster and stabler ways to compute
% this (balance and then exponentiate...
for t = T
    expmXp = expm( X * t );
    helper6 = expmXp * xvec6;
    mesaout6 = zeros( NSamplesPlot );
    for n = 1:NSamplesPlot
        for l = 1:NSamplesPlot
            mesaout6( n, l ) = (4*pi^2)^-1 * sum( helper6 .* ( exp( 1i * u( n, l ) * multi_idx( :, 1 ) ) .* exp( 1i * v( n, l ) * multi_idx( :, 2 ) ) ) );
        end
    end
    h = figure; hold on;
    contour( u, v, mesaout6 .* conj( mesaout6 ) ); colorbar;
    title(t);
    xlabel('x1');
    ylabel('x2');
    shading interp;
    plot( Y( :, 1 ), Y( :, 2 ), 'k', 'LineWidth', 2 );
%     filename = sprintf( 'DGHD_evolution/evolution_%3.25f.pdf', t );
%     print( h, '-dpdf', filename );
end
