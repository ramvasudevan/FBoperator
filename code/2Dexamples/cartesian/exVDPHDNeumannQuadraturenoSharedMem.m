%% parameters
NBases = 51;
compute_nulls = 0;
NSamplesInt = 125;
NSamplesPlot = 101;
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
[ u, v ] = meshgrid( linspace( -pi, pi, NSamplesPlot ) );
u = u( : );
v = v( : );
%% compute functions for later use
tic
[ f1, f2, f3, f4 ] = vdpR2( tx, ty );
PlotTrig = zeros( size( multi_idx, 1 ), NSamplesPlot^2 );
Trig = zeros( size( multi_idx, 1 ), NSamplesInt^2 );
DfTrig = zeros( size( multi_idx, 1 ), NSamplesInt^2 );
nc = ones( size( multi_idx, 1 ), 1 )/sqrt( pi );
nc( 1 ) = 1/sqrt( 2 * pi );
pc = ones( size( multi_idx, 1 ), 1 )/sqrt( pi );
pc( 1 ) = 1/sqrt( 2 * pi );
for k = 1:size( multi_idx, 1 )
    cx = multi_idx( k, 1 );
    cy = multi_idx( k, 2 );
    if( mod( cx, 2 ) && mod( cy, 2 ) )
        PlotTrig( k, : ) = sin( cx * u/2 ) .* sin( cy * v/2 ) .* pc( cx + 1 ) .* pc( cy + 1 );
        Trig( k, : ) = sin( cx * tx/2 ) .* sin( cy * ty/2 ) .* nc( cx + 1 ) .* nc( cy + 1 );
        DfTrig( k, : ) =  Weights .* ( f1 .* cx/2 .* cos( cx * tx/2 ) .* sin( cy * ty/2 ) + f2 .* cy/2 .* sin( cx * tx/2 ) .* cos( cy * ty/2 ) + ...
            1/2 .* ( f3 + f4 ) .* sin( cx * tx/2 ) .* sin( cy * ty/2 ) ) .* nc( cx + 1 ) .* nc( cy + 1 );
    elseif( mod( cx, 2 ) )
        PlotTrig( k, : ) = sin( cx * u/2 ) .* cos( cy * v/2 ) .* pc( cx + 1 ) .* pc( cy + 1 );
        Trig( k, : ) = sin( cx * tx/2 ) .* cos( cy * ty/2 ) .* nc( cx + 1 ) .* nc( cy + 1 );
        DfTrig( k, : ) =  Weights .* ( f1 .* cx/2 .* cos( cx * tx/2 ) .* cos( cy * ty/2 ) + f2 .* -cy/2 .* sin( cx * tx/2 ) .* sin( cy * ty/2 ) + ...
            1/2 .* ( f3 + f4 ) .* sin( cx * tx/2 ) .* cos( cy * ty/2 ) ) .* nc( cx + 1 ) .* nc( cy + 1 );
    elseif( mod( cy, 2 ) )
        PlotTrig( k, : ) = cos( cx * u/2 ) .* sin( cy * v/2 ) .* pc( cx + 1 ) .* pc( cy + 1 );
        Trig( k, : ) = cos( cx * tx/2 ) .* sin( cy * ty/2 ) .* nc( cx + 1 ) .* nc( cy + 1 );
        DfTrig( k, : ) = Weights .* ( f1 .* -cx/2 .* sin( cx * tx/2 ) .* sin( cy * ty/2 ) + f2 .* cy/2 .* cos( cx * tx/2 ) .* cos( cy * ty/2 ) + ...
            1/2 .* ( f3 + f4 ) .* cos( cx * tx/2 ) .* sin( cy * ty/2 ) ) .* nc( cx + 1 ) .* nc( cy + 1 );
    else
        PlotTrig( k, : ) = cos( cx * u/2 ) .* cos( cy * v/2 ) .* pc( cx + 1 ) .* pc( cy + 1 );
        Trig( k, : ) = cos( cx * tx/2 ) .* cos( cy * ty/2 ) .* nc( cx + 1 ) .* nc( cy + 1 );
        DfTrig( k, : ) = Weights .* ( f1 .* -cx/2 .* sin( cx * tx/2 ) .* cos( cy * ty/2 ) + f2 .* -cy/2 .* cos( cx * tx/2 ) .* sin( cy * ty/2 ) + ...
            1/2 .* ( f3 + f4 ) .* cos( cx * tx/2 ) .* cos( cy * ty/2 ) ) .* nc( cx + 1 ) .* nc( cy + 1 );
    end
end
toc
%% compute the eigenrepresentation of X
% in the complex exponential basis the indices correspond to the power of
% the exponential, for the trigonometric basis the even ones correspond to
% cos and the odd ones correspond to the sin one.
% tic
% X = zeros( size( multi_idx, 1 ) );
% for k = 1:size( multi_idx, 1 )
%     X( k, : ) = sum( bsxfun( @times, squeeze( DfTrig( k, : ) ), Trig ), 2 )';
% end
% toc
tic
X = DfTrig * Trig';
toc
X( abs( X ) < 1e-10 ) = 0;
X = sparse( X );

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
    %         toc
    %             tic
    %             [ U, D, V ] = svd( X );
    %             zeroeig = find( abs( real( diag( D ) ) ) < 1e-10 )
    %             toc
    tic
    zeroeig = 1:10;
    [ V, D ] = eigs( X, length( zeroeig ), 'SM' );
    eigens = diag( D )
    toc;
end

%% just for plotting
tic
out = zeros( length( zeroeig ), NSamplesPlot, NSamplesPlot );
for m = 2:2:length( zeroeig )
    out( m, :, : ) = Neumann2DPlot( V( :, zeroeig( m ) ), NSamplesPlot, PlotTrig );
end
toc

%% plot vector field
[ Du, Dv, ~, ~ ] = vdpR2( u, v );
[ TY, Y ] = ode23( @vdpR2_plot, [0 0.75], [ -2.5 -2.5 ] );
% figure; hold on;
% quiver( u, v, Du, Dv );
% plot( Y( :, 1 ), Y( :, 2 ), 'g' );
% print( h, '-dpdf', 'indicator_mesa_images/vf_and_lc.pdf' );

%% actual plotting
u = reshape( u, [ NSamplesPlot NSamplesPlot ] );
v = reshape( v, [ NSamplesPlot NSamplesPlot ] );
for m = 2:2:length( zeroeig )
    h = figure; hold on;
    %     set(gcf,'renderer','painters')
    surf( u, v, squeeze( out( m, :, : ) ).* conj( squeeze( out( m, :, : ) ) ) ); colorbar;
    %     contour( u, v, squeeze( out( m, :, : ) ).* conj( squeeze( out( m, :, : ) ) ), 100 ); colorbar;
    shading flat
    plot( Y( :, 1 ), Y( :, 2 ), 'g', 'LineWidth', 10 );
    title(imag(eigens(zeroeig(m))));
    xlabel('x1');
    ylabel('x2');
    %     filename = sprintf( 'vdp_densities/eigdist_image_%3.10f.pdf', imag( eigens( m ) ) );
    %     print( h, '-dpdf', filename );
end

%% plot mesa function evolution
syms x1 x2;
[ TY, Y ] = ode23( @vdpR2_plot, [0 0.75], [ 2 -2 ] );
H = matlabFunction( ( heaviside( x1 - 1.75 ) - heaviside( x1 - 2.25 ) ) * ( heaviside( x2 + 2.25 ) - heaviside( x2 + 1.75 ) ) );
% H = matlabFunction( ( heaviside( x1 + 0.5 ) - heaviside( x1 - 0.5 ) ) * ( heaviside( x2 + 1/2 ) - heaviside( x2 - 1/2 ) ) );
F = Weights .* H( tx, ty );
evo0 = Trig * ( Weights .* H( tx, ty ) );

% length of time to plot
T = linspace( 0, 0.75, 30 );

% consider doing the following: [ V, D ] = eig( X ); and expm( X * t ) = V * diag( exp( diag( t * D ) ) ) * V'
for t = T
    expmXp = expm( X * t );
    evo = expmXp * evo0;
    evoPlot = Neumann2DPlot( evo, NSamplesPlot, PlotTrig );
    
    figure; hold on;
    surf( u, v, evoPlot .* conj( evoPlot ) ); colorbar;
    title( t );
    xlabel('x1');
    ylabel('x2');
    shading flat;
%     view( 0, 89 );
    plot( Y( TY < t, 1 ), Y( TY < t, 2 ), 'g', 'LineWidth', 10 );
    set(gcf,'Renderer','Zbuffer')
end
