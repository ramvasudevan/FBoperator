%% parameters
NBases = 5;
compute_all_eigs = 1;

%% generate multi-idx to make clear what each element in row is
multi_idx = zeros( ( 2 * NBases + 1 )^3, 3 );
counter = 0;
for n = -NBases:NBases
    for m = -NBases:NBases
        for l = -NBases:NBases
            counter = counter + 1;
            multi_idx( counter, : ) = [ n m l ];
        end
    end
end

%% compute p and its DFT
Fs = 1000;
t = 0:1/Fs:( 2 * pi );
p = t;
nfft = length( p );
% Take fft, padding with zeros so that length(X) is equal to nfft
P = fft( p, nfft );
P = fftshift( P );
% Frequency vector
freqs = (-nfft/2:nfft/2-1)*Fs/nfft;
% p's representation in our bases
pvec = zeros( 2 * NBases + 1, 1 );
for k = -NBases:NBases
    search_freq = k/( 2 * pi );
    if ( find( freqs == search_freq ) )
        pvec( k + NBases + 1 ) = P( freqs == search_freq );
    else
        min_idx = find( freqs < search_freq , 1, 'last' );
        max_idx = find( freqs > search_freq , 1 );
        min_freq = freqs( min_idx );
        max_freq = freqs( max_idx );
        t_eval = ( search_freq - min_freq )/( max_freq - min_freq );
        pvec( k + NBases + 1 ) = P( max_idx ) * t_eval + P( min_idx ) * ( 1 - t_eval );
    end
end

%% compute f(\theta) and its DFT
Fs = 1000;
t = ( -pi ):1/Fs:( pi );
f = zeros( size( t ) );
f( t < pi/2 ) = exp( ( t( t < pi/2 ).^2 - pi^2/4 ).^( -1 ) );
f( t < -pi/2 ) = 0;
nfft = length( f );
% Take fft, padding with zeros so that length(X) is equal to nfft
F = fft( f, nfft );
F = fftshift( F );
% Frequency vector
freqs = (-nfft/2:nfft/2-1)*Fs/nfft;
% f's representation in our bases
fvec = zeros( 2 * NBases + 1, 1 );
for k = -NBases:NBases
    search_freq = k/( 2 * pi );
    if ( find( freqs == search_freq ) )
        fvec( k + NBases + 1 ) = F( freqs == search_freq );
    else
        min_idx = find( freqs < search_freq , 1, 'last' );
        max_idx = find( freqs > search_freq , 1 );
        min_freq = freqs( min_idx );
        max_freq = freqs( max_idx );
        t_eval = ( search_freq - min_freq )/( max_freq - min_freq );
        fvec( k + NBases + 1 ) = F( max_idx ) * t_eval + F( min_idx ) * ( 1 - t_eval );
    end
end

%% compute the eigenrepresentation of X
% t portion of representation
Xt = zeros( ( 2 * NBases + 1 )^3, ( 2 * NBases + 1 )^3  );
for k = 1:size( multi_idx, 1 )
    Xt( k, k ) = -1i * multi_idx( k, 1 );
end
% Xt = sparse( 1:( ( 2 * NBases + 1 )^3 ), 1:( ( 2 * NBases + 1 )^3 ), double( multi_idx( :, 1 ) ) );

% theta portion of representation
Xtheta = zeros( ( 2 * NBases + 1 )^3, ( 2 * NBases + 1 )^3 );
for k = 1:size( multi_idx, 1 )
    c_idx = multi_idx( k, : );
    for j = -NBases:NBases
        h_idx = ismember( multi_idx, [ c_idx( 1:2 ) c_idx( 3 ) + j ], 'rows' );
        Xtheta( h_idx, k ) = -1i * c_idx( 2 ) * pvec( j + NBases + 1 );
    end
end

% p portion of representation
Xp = zeros( ( 2 * NBases + 1 )^3, ( 2 * NBases + 1 )^3 );
for k = 1:size( multi_idx, 1 )
    c_idx = multi_idx( k, : );
    
    %     % t-portion
    %     for j = -NBases:NBases
    %         h_idx = ismember( multi_idx, [ c_idx( 1 ) - 2 c_idx( 2 ) + j c_idx( 3 ) ], 'rows' );
    %         Xp( h_idx, k ) = Xp( h_idx, k ) - c_idx( 3 )/2 * fvec( j + NBases + 1 );
    %         h_idx = ismember( multi_idx, [ c_idx( 1 ) + 2 c_idx( 2 ) + j c_idx( 3 ) ], 'rows' ) ;
    %         Xp( h_idx, k ) = Xp( h_idx, k ) + c_idx( 3 )/2 * fvec( j + NBases + 1 );
    %     end
    
    % theta-portion
    %     h_idx = ismember( multi_idx, [ c_idx( 1 ) c_idx( 2 ) - 2 c_idx( 3 ) ], 'rows' );
    %     Xp( h_idx, k ) = -c_idx( 3 )/2;
    %     h_idx = ismember( multi_idx, [ c_idx( 1 ) c_idx( 2 ) + 2 c_idx( 3 ) ], 'rows' );
    %     Xp( h_idx, k ) = c_idx( 3 )/2;
    h_idx = ismember( multi_idx, [ c_idx( 1 ) c_idx( 2 ) - 2 c_idx( 3 ) ], 'rows' );
    Xp( h_idx, k ) = -c_idx( 3 )/2;
    h_idx = ismember( multi_idx, [ c_idx( 1 ) c_idx( 2 ) + 2 c_idx( 3 ) ], 'rows' );
    Xp( h_idx, k ) = c_idx( 3 )/2;
    
    % p-portion
    %     Xp( k, k ) = ( 1 + 1i * c_idx( 3 ) * pvec( c_idx( 3 ) + NBases + 1 ) );
    %     Xp( k, k ) = Xp( k, k ) - c_idx( 2 ) * c_idx( 3 );
end

%% compute eigenvalues of different portions
if ( compute_all_eigs )
    %     [ Vt, Dt ] = eig( Xt, eye( size( Xt ) ) );
    %     [ Vp, Dp ] = eig( Xp, eye( size( Xp ) ) );
    %     [ Vtheta, Dtheta ] = eig( Xtheta, eye( size( Xtheta ) ) );
    [ ~, helpert ] = spspaces( Xt, 2 );
    Vt = helpert{1};
    zeroeigt = helpert{3};
    [ ~, helperp ] = spspaces( Xp, 2 );
    Vp = helperp{1};
    zeroeigp = helperp{3};
    [ ~, helpertheta ] = spspaces( Xtheta, 2 );
    Vttheta = helpertheta{1};
    zeroeigtheta = helpertheta{3};
    
else
    % this works fast
    [ Vt, Dt ] = eigs( Xt, eye( size( Xt ) ), 1, 1e-10 );
    [ Vp, Dp ] = eigs( Xp, eye( size( Xp ) ), 1, 1e-10 );
    [ Vtheta, Dtheta ] = eigs( Xtheta, eye( size( Xtheta ) ), 1, 1e-10 );
    
    zeroeigt = 1;
    zeroeigp = 1;
    zeroeigtheta = 1;
end


%% just for plotting
NSamples = 1001;
x = linspace( 0, 2 * pi, NSamples );
outt = zeros( 1, NSamples );
outp = zeros( 1, NSamples );
outtheta = zeros( 1, NSamples );
for j = 1:NSamples
    for k = 1:size( multi_idx, 1 )
        outt( j ) = outt( j ) + Vt( k, zeroeigt( 1 ) ) * exp( 1i * multi_idx( k, 1 ) * x( j ) ) ...
            * exp( 1i *  multi_idx( k, 2 ) * x( j ) ) * exp( 1i * multi_idx( k, 3 ) * x( j ) );
        outp( j ) = outp( j ) + Vp( k, zeroeigp( 1 ) ) * exp( 1i * multi_idx( k, 1 ) * x( j ) ) ...
            * exp( 1i *  multi_idx( k, 2 ) * x( j ) ) * exp( 1i * multi_idx( k, 3 ) * x( j ) );
        outtheta( j ) = outtheta( j ) + Vtheta( k, zeroeigtheta( 1 ) ) * exp( 1i * multi_idx( k, 1 ) * x( j ) ) ...
            * exp( 1i *  multi_idx( k, 2 ) * x( j ) ) * exp( 1i * multi_idx( k, 3 ) * x( j ) );
    end
end


figure;
subplot( 3, 1, 1 ); plot( x, outt );
subplot( 3, 1, 2 ); plot( x, outtheta );
subplot( 3, 1, 3 ); plot( x, outp );


