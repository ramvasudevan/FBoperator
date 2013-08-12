%% parameters
NBases = 1;
compute_nulls = 1;

%% generate multi-idx to make clear what each element in row is
[ gridx, gridy ] = meshgrid( -NBases:NBases );
multi_idx = zeros( size( gridx, 1 ) * size( gridx, 2 ), 2 );
multi_idx( :, 1 ) = gridx( : );
multi_idx( :, 2 ) = gridy( : );

%% compute f_{1,2} and its DFT 
%constructing the signals
Fs = 512;
t = 0:( ( 2 * pi )/Fs ):( ( 2 * pi ) - ( 2 * pi)/Fs );
% t = ( -2 * pi ):( ( 4 * pi )/Fs ):( ( 2 * pi ) - ( 4 * pi)/Fs );
f1 = zeros( length( t ) );
f2 = zeros( length( t ) );
f3 = zeros( length( t ) );
f4 = zeros( length( t ) );
for j = 1:length( t )
    for k = 1:length( t )
        f1( k, j ) = pi * cos( t( j ) ) * cos( t( k ) );
        f2( k, j ) = pi * sin( t( j ) ) * cos( t( k ) );
        f3( k, j ) = pi * cos( t( j ) ) * cos( t( k ) );
        f4( k, j ) = pi * cos( t( j ) ) * sin( t( k ) );
    end
end
% computing their FFTs
nfft = length( t );
F1 = fftshift( fft2( f1, nfft, nfft ) )/nfft^2;
F2 = fftshift( fft2( f2, nfft, nfft ) )/nfft^2;
F3 = fftshift( fft2( f3, nfft, nfft ) )/nfft^2;
F4 = fftshift( fft2( f4, nfft, nfft ) )/nfft^2;
% frequency vector
freqs = ( -nfft/2:( nfft/2 - 1 ) ) * Fs/nfft;
sampled_freqs = -NBases:NBases;
[ freqsx, freqsy ] = meshgrid( freqs );
f1vec = interp2( freqsx, freqsy, F1, gridx, gridy, 'cubic' );
f2vec = interp2( freqsx, freqsy, F2, gridx, gridy, 'cubic' );
f3vec = interp2( freqsx, freqsy, F3, gridx, gridy, 'cubic' );
f4vec = interp2( freqsx, freqsy, F4, gridx, gridy, 'cubic' );

%% compute the eigenrepresentation of X
tic
% X1 = zeros( size( multi_idx, 1 ) );
% X2 = zeros( size( multi_idx, 1 ) );
X = zeros( size( multi_idx, 1 ) );
for k = 1:size( multi_idx, 1 )
    c_idx = multi_idx( k, : );
    s_idx( :, 1 ) = c_idx( 1 ) + gridx( : );
    s_idx( :, 2 ) = c_idx( 2 ) + gridy( : );
    f_idx = intersect( intersect( find( c_idx( 1 ) + gridx( : ) >= -NBases ), find( c_idx( 1 ) + gridx( : ) <= NBases ) ), ...
        intersect( find( c_idx( 2 ) + gridy( : ) >= -NBases ), find( c_idx( 2 ) + gridy( : )  <= NBases ) ) );
    h_idx = any( all( bsxfun( @eq, reshape( s_idx.', 1, 2, [] ), multi_idx), 2 ), 3 );
%     X1( h_idx, k ) = f1vec( f_idx ) + 1i * c_idx( 1 ) * f2vec( f_idx );
%     X2( h_idx, k ) = -f3vec( f_idx ) - 1i * c_idx( 2 ) * f4vec( f_idx );
    X( h_idx, k ) = f1vec( f_idx ) - f3vec( f_idx ) + ...
        1i * ( c_idx( 1 ) * f2vec( f_idx ) - c_idx( 2 ) * f4vec( f_idx ) );
end
toc

%% zero out small entries
X( abs( X ) < 1e-10 ) = 0;

%% compute eigenvalues of different portions
if ( compute_nulls )
    tic
    [ ~, helper ] = spspaces( X, 2 );
    [ ~, helper2 ] = spspaces( X2, 2 );
    toc
    V = helper{1};
    zeroeig = helper{3};
else
    tic
    [ V, D ] = eig( X, eye( size( X ) ) );
    toc
    zeroeig = [];
    for k = 1:length( D )
        if ( abs( D( k, k ) ) == 0 )
            zeroeig = [ zeroeig; k ];
        end
    end
end

%% generate the LP to find the distribution using Parseval's Theorem ( first set correspond to positivity requirements,
% last set correspond to real-valuedness requirement)
tic
Aineq = zeros( 2 * ( NBases + 1 )^2 + 1 - 2 + 2 * NBases^2, length( zeroeig ) );
Aeq = zeros( ( 2 * NBases + 1 )^2, length( zeroeig ) ); % last constraint ensures its a probability distribution and LOTS OF additional constraints here but whatever
bineq = zeros( 2 * ( NBases + 1 )^2 + 1 - 2 + 2 * NBases^2, 1 );
beq = zeros( ( 2 * NBases + 1 )^2, 1 ); % lots of additional constraints here but whatever.
beq( end ) = 1/( 2 * pi ); % last constraint ensures its a probability distribution
obj = zeros( length( zeroeig ), 1 ); % just check feasibility

z_idx = any( all( bsxfun( @eq, reshape( [ 0 0 ].', 1, 2, [] ), multi_idx), 2 ), 3 );
% make the square integral of the function equal to the integral of the
% function
H = diag( sum( V( :, zeroeig ) .* V( :, zeroeig ) ) );
obj( : ) = -1/2 * V( z_idx, zeroeig );

% take care of the constant term
Aineq( 1, : ) = V( z_idx, zeroeig );

% all possible combinations
comb_idx01 = [ 1 1; 1 -1 ];
sine_idx01 = repmat( [ -1; 1 ], [ 1 length( zeroeig ) ] );
comb_idx10 = [ 1 1; -1 1 ];
sine_idx10 = repmat( [ -1; 1 ], [ 1 length( zeroeig ) ] );
comb_idx11 = [ 1 1; 1 -1; -1 1; -1 -1 ];
sine_idx11 = repmat( [ 1; -1; -1; 1 ], [ 1 length( zeroeig ) ] );

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
            Aineq( r_counter, : ) = V( z_idx, zeroeig ) + 0.5 * sum( V( h_idx, zeroeig ) );
        else
            Aineq( r_counter, : ) = V( z_idx, zeroeig ) + 0.25 * sum( V( h_idx, zeroeig ) );
        end
        r_counter = r_counter + 1;
        
        % the 1 + sine terms
        if ( j == 0 || k == 0 )
            Aineq( r_counter, : ) = V( z_idx, zeroeig ) - 0.5 * sum( sine_idx .* V( h_idx, zeroeig ) );
        else
            Aineq( r_counter, : ) = V( z_idx, zeroeig ) - 0.25 * sum( sine_idx .* V( h_idx, zeroeig ) );
        end
        r_counter = r_counter + 1;
    end
end

% all possible combinations
comb_idx = [ 1 1; 1 -1; -1 1; -1 -1 ];
sc_idx = repmat( [ -1; -1; 1; 1 ], [ 1 length( zeroeig ) ] );
cs_idx = repmat( [ -1; 1; -1; 1 ], [ 1 length( zeroeig ) ] );

for k = 1:NBases
    for j = 1:NBases
        s_idx = repmat( [ j k ], [ 4 1 ] ) .* comb_idx;
        h_idx = zeros( size( s_idx, 1 ), 1 );
        for l = 1:length( h_idx )
            h_idx( l ) = find( any( all( bsxfun( @eq, reshape( s_idx( l, : ).', 1, 2, [] ), multi_idx), 2 ), 3 ) );
        end
        
        % the 1 + sin(x) * cos(y)
        Aineq( r_counter, : ) = V( z_idx, zeroeig ) + 0.25 * 1i * sum( sc_idx .* V( h_idx, zeroeig ) );
        r_counter = r_counter + 1;
        
        % the 1 + cos(x) * sin(y)
        Aineq( r_counter, : ) = V( z_idx, zeroeig ) + 0.25 * 1i * sum( cs_idx .* V( h_idx, zeroeig ) );
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
        
        Aeq( r_counter, : ) = sum( repmat( [ 1; -1 ], [ 1 length( zeroeig ) ] ).* V( h_idx, zeroeig ) );
        r_counter = r_counter + 1;
    end
end

% ensure its a probability density
Aeq( r_counter, : ) = V( z_idx, zeroeig ) * ( 2 * pi );
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
            out( n, l ) = out( n, l ) +  sum( ( scale .* V( :, zeroeig ) )' * helper( :, l, n ) );
        end
    end
    toc
    figure; mesh( x, x, real( out ) ); colorbar;
end
