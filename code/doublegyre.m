%% parameters
NBases = 1;
compute_nulls = 1;

%% generate multi-idx to make clear what each element in row is
[ gridx, gridy ] = meshgrid( -NBases:NBases );
multi_idx = zeros( size( gridx, 1 ) * size( gridx, 2 ), 2 );
multi_idx( :, 1 ) = gridx( : );
multi_idx( :, 2 ) = gridy( : );

%% compute the eigenrepresentation of X
tic
X1 = zeros( size( multi_idx, 1 ) );
X2 = zeros( size( multi_idx, 1 ) );
f1vec = zeros( size( multi_idx, 1 ), 1 );
f2vec = zeros( size( multi_idx, 1 ), 1 );
s_idx = zeros( 4, 2 );
for k = 1:size( multi_idx, 1 )
    c_idx = multi_idx( k, : );
    s_idx( :, 1 ) = c_idx( 1 ) +  [ 1 1 -1 -1 ];
    s_idx( :, 2 ) = c_idx( 2 ) + [ 1 -1 1 -1 ];
    h_idx11 = any( all( bsxfun( @eq, reshape( s_idx( 1, : ).', 1, 2, [] ), multi_idx), 2 ), 3 );
    h_idx1m1 = any( all( bsxfun( @eq, reshape( s_idx( 2, : ).', 1, 2, [] ), multi_idx), 2 ), 3 );
    h_idxm11 = any( all( bsxfun( @eq, reshape( s_idx( 3, : ).', 1, 2, [] ), multi_idx), 2 ), 3 );
    h_idxm1m1 = any( all( bsxfun( @eq, reshape( s_idx( 4, : ).', 1, 2, [] ), multi_idx), 2 ), 3 );
    if ( find( h_idx11 ) )
        f1vec( h_idx11 ) = pi/4 * ( 1 + c_idx( 1 ) );
        f2vec( h_idx11 ) = pi/4 * ( -1 - c_idx( 2 ) );
    end
    if ( find( h_idx1m1 ) )
        f1vec( h_idx1m1 ) = pi/4 * ( 1 + c_idx( 1 ) );
        f2vec( h_idx1m1 ) = pi/4 * ( -1 + c_idx( 2 ) );
    end
    if ( find( h_idxm11 ) )
        f1vec( h_idxm11 ) = pi/4 * ( 1 - c_idx( 1 ) );
        f2vec( h_idxm11 ) = pi/4 * ( -1 - c_idx( 2 ) );
    end
    if ( find( h_idxm1m1 ) )
        f1vec( h_idxm1m1 ) = pi/4 * ( 1 - c_idx( 1 ) );
        f2vec( h_idxm1m1 ) = pi/4 * ( -1 + c_idx( 2 ) );
    end
    h_idx = any( all( bsxfun( @eq, reshape( s_idx.', 1, 2, [] ), multi_idx), 2 ), 3 );
    X1( h_idx, k ) = f1vec( h_idx ) +  f2vec( h_idx );
%     X1( h_idx, k ) = f1vec( h_idx );
%     X2( h_idx, k ) = f2vec( h_idx );
end
toc

%% zero out small entries
% X1( abs( X1 ) < 1e-10 ) = 0;
% X2( abs( X2 ) < 1e-10 ) = 0;

%% compute eigenvalues of different portions
if ( compute_nulls )
    tic
    [ ~, helper1 ] = spspaces( X1, 2 );
    [ ~, helper2 ] = spspaces( X2, 2 );
    toc
    V1 = helper1{1};
    zeroeig1 = helper1{3};
    V2 = helper2{1};
    zeroeig2 = helper2{3};
else
    % DO NOT USE ITS AWFULLY INACCURATE!
    tic
    [ V1, D1 ] = eig( X1, eye( size( X1 ) ) );
    [ V2, D2 ] = eig( X2, eye( size( X2 ) ) );
    toc
    zeroeig1 = [];
    zeroeig2 = [];
    for k = 1:length( D1 )
        if ( abs( D1( k, k ) ) == 0 )
            zeroeig1 = [ zeroeig1; k ];
        end
        if ( abs( D2( k, k ) ) == 0 )
            zeroeig2 = [ zeroeig2; k ];
        end
    end
end

% %% just for plotting
% NSamples = 101;
% tic
% x = linspace( -pi, pi, NSamples );
% [ u, v ] = meshgrid( x );
% helper = reshape( exp( 1i * kron( u( : ), multi_idx( :, 1 ) ) ).* exp( 1i * kron( v( : ), multi_idx( :, 2 ) ) ), ...
%     size( multi_idx, 1 ), NSamples, NSamples );
% out2 = zeros( length( zeroeig2 ), NSamples, NSamples );
% parfor m = 1:length( zeroeig1 )
%     for n = 1:NSamples
%         for l = 1:NSamples
%             out1( m, n, l ) = V1( :, zeroeig1( m ) )' * helper( :, l, n );
%         end
%     end
% end
% toc
% %% actual plotting
% for m = 1:size( zeroeig1, 1 )
% % for m = 75
%     figure;
%     mesh( x, x, real( squeeze( out1( m, :, : ) ) ) ); colorbar;
%     pause;
% end
%
% % figure; hold on;
% for m = 1:size( zeroeig2, 1 )
% % for m = 75
%     figure;
%     mesh( x, x, real( squeeze( out2( m, : , : ) ) ) ); colorbar;
%     pause;
% end


%% generate the LP to find the distribution using Parseval's Theorem ( first set correspond to positivity requirements,
% last set correspond to real-valuedness requirement)
tic
Aineq = zeros( 2 * ( NBases + 1 )^2 + 1 - 2 + 2 * NBases^2, length( zeroeig1 ) );
Aeq = zeros( ( 2 * NBases + 1 )^2, length( zeroeig1 ) ); % last constraint ensures its a probability distribution and LOTS OF additional constraints here but whatever
bineq = zeros( 2 * ( NBases + 1 )^2 + 1 - 2 + 2 * NBases^2, 1 );
beq = zeros( ( 2 * NBases + 1 )^2, 1 ); % lots of additional constraints here but whatever.
beq( end ) = 1/( 2 * pi ); % last constraint ensures its a probability distribution
obj = zeros( length( zeroeig1 ), 1 ); % just check feasibility

z_idx = any( all( bsxfun( @eq, reshape( [ 0 0 ].', 1, 2, [] ), multi_idx), 2 ), 3 );
% make the square integral of the function equal to the integral of the
% function
H = diag( sum( V1( :, zeroeig1 ) .* V1( :, zeroeig1 ) ) );
obj( : ) = -1/2 * V1( z_idx, zeroeig1 );

% take care of the constant term
Aineq( 1, : ) = V1( z_idx, zeroeig1 );

% all possible combinations
comb_idx01 = [ 1 1; 1 -1 ];
sine_idx01 = repmat( [ -1; 1 ], [ 1 length( zeroeig1 ) ] );
comb_idx10 = [ 1 1; -1 1 ];
sine_idx10 = repmat( [ -1; 1 ], [ 1 length( zeroeig1 ) ] );
comb_idx11 = [ 1 1; 1 -1; -1 1; -1 -1 ];
sine_idx11 = repmat( [ 1; -1; -1; 1 ], [ 1 length( zeroeig1 ) ] );

% ensure positivity
r_counter = 2;
for j = 0:NBases
    for k = 0:NBases
        if ( j == 0 && k == 0 )
            continue;
        end
        
        if ( j == 0 )
            s_idx = repmat( [ j k ], [ 2 1 ] ) .* comb_idx01;
            sine_idx = sine_idx01;
        elseif( k == 0 )
            s_idx = repmat( [ j k ], [ 2 1 ] ) .* comb_idx10;
            sine_idx = sine_idx10;
        else
            s_idx = repmat( [ j k ], [ 4 1 ] ) .* comb_idx11;
            sine_idx = sine_idx11;
        end
        
        %         s_idx = repmat( [ j k ], [ 4, 1 ] ) .* comb_idx;
        h_idx = zeros( size( s_idx, 1 ), 1 );
        for l = 1:length( h_idx )
            h_idx( l ) = find( any( all( bsxfun( @eq, reshape( s_idx( l, : ).', 1, 2, [] ), multi_idx), 2 ), 3 ) );
        end
        
        % the 1 + cosine terms
        if ( j == 0 || k == 0 )
            Aineq( r_counter, : ) = V1( z_idx, zeroeig1 ) + 0.5 * sum( V1( h_idx, zeroeig1 ) );
        else
            Aineq( r_counter, : ) = V1( z_idx, zeroeig1 ) + 0.25 * sum( V1( h_idx, zeroeig1 ) );
        end
        r_counter = r_counter + 1;
        
        % the 1 + sine terms
        if ( j == 0 || k == 0 )
            Aineq( r_counter, : ) = V1( z_idx, zeroeig1 ) - 0.5 * sum( sine_idx .* V1( h_idx, zeroeig1 ) );
        else
            Aineq( r_counter, : ) = V1( z_idx, zeroeig1 ) - 0.25 * sum( sine_idx .* V1( h_idx, zeroeig1 ) );
        end
        r_counter = r_counter + 1;
    end
end

% all possible combinations
comb_idx = [ 1 1; 1 -1; -1 1; -1 -1 ];
sc_idx = repmat( [ -1; -1; 1; 1 ], [ 1 length( zeroeig1 ) ] );
cs_idx = repmat( [ -1; 1; -1; 1 ], [ 1 length( zeroeig1 ) ] );

for k = 1:NBases
    for j = 1:NBases
        s_idx = repmat( [ j k ], [ 4 1 ] ) .* comb_idx;
        h_idx = zeros( size( s_idx, 1 ), 1 );
        for l = 1:length( h_idx )
            h_idx( l ) = find( any( all( bsxfun( @eq, reshape( s_idx( l, : ).', 1, 2, [] ), multi_idx), 2 ), 3 ) );
        end
        
        % the 1 + sin(x) * cos(y)
        Aineq( r_counter, : ) = V1( z_idx, zeroeig1 ) + 0.25 * 1i * sum( sc_idx .* V1( h_idx, zeroeig1 ) );
        r_counter = r_counter + 1;
        
        % the 1 + cos(x) * sin(y)
        Aineq( r_counter, : ) = V1( z_idx, zeroeig1 ) + 0.25 * 1i * sum( cs_idx .* V1( h_idx, zeroeig1 ) );
        r_counter = r_counter + 1;
    end
end

% ensure real-valuedness
r_counter = 1;
for j = -NBases:NBases
    for k = -NBases:NBases
        if ( j == 0 && k == 0 )
            continue;
        end
        h_idx = any( all( bsxfun( @eq, reshape( [ j k; -j -k ].', 1, 2, [] ), multi_idx), 2 ), 3 );
        
        Aeq( r_counter, : ) = sum( repmat( [ 1; -1 ], [ 1 length( zeroeig1 ) ] ).* V1( h_idx, zeroeig1 ) );
        r_counter = r_counter + 1;
    end
end

% ensure its a probability density
Aeq( r_counter, : ) = V1( z_idx, zeroeig1 ) * ( 2 * pi );
toc
tic;
% [ result, value, exitflag, output]  = cplexlp( obj, -Aineq, -bineq, Aeq, beq );
[ result, value, exitflag, output]  = cplexqp( H, obj, -Aineq, -bineq, Aeq, beq );
toc

%% plot the element in the nullspace that comes out of the optimization
if ( exitflag == -2 )
    disp(' No feasible solution' );
elseif( exitflag == -8 )
    disp( 'Unbounded solution' );
elseif( exitflag == 1 )
    NSamples = 201;
    tic
    x = linspace( -pi, pi, NSamples );
    [ u, v ] = meshgrid( x );
    scale = repmat( result', [ size( multi_idx, 1 ), 1 ] );
    helper = reshape( exp( 1i * kron( u( : ), multi_idx( :, 1 ) ) ).* exp( 1i * kron( v( : ), multi_idx( :, 2 ) ) ), ...
        size( multi_idx, 1 ), NSamples, NSamples );
    out = zeros( NSamples, NSamples );
    parfor n = 1:NSamples
        for l = 1:NSamples
            out( n, l ) = out( n, l ) +  sum( ( scale .* V1( :, zeroeig1 ) )' * helper( :, l, n ) );
        end
    end
    toc
    figure; mesh( x, x, real( out ) ); colorbar;
end

% %% plot vector field
% [ u, v ] = meshgrid( -pi:0.1*pi:( pi ), -pi:0.1*pi:( pi ) );
% Du = zeros( size( u ) );
% Dv = zeros( size( v ) );
% for j = 1:size( u, 1 )
%     for k = 1:size( u, 2 )
%         Du( j, k ) = -pi * sin( u( j, k ) ) * cos( v( j, k ) );
%         Dv( j, k ) = pi * cos( u( j, k ) ) * sin( v( j, k ) );
%     end
% end
% figure; quiver( u, v, Du, Dv )


