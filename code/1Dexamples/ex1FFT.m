%% parameters
NBases = 100;
compute_nulls = 1;

%% generate multi-idx to make clear what each element in row is
multi_idx = (-NBases:NBases)';

%% compute f_{1,2} and its DFT
% constructing the signals
Fs = 5000;
t = 0:( ( 2 * pi )/Fs ):( ( 2 * pi ) - ( 2 * pi)/Fs );
f1 = cos( t );
f2 = sin( t );
% computing their FFTs
nfft = length( t );
F1 = fftshift( fft( f1, nfft ) )/nfft;
F2 = fftshift( fft( f2, nfft ) )/nfft;
% frequency vector
freqs = ( -nfft/2:( nfft/2 - 1 ) ) * Fs/nfft;
sampled_freqs = -NBases:NBases;
f1vec = interp1( freqs, F1, sampled_freqs, 'cubic' );
f2vec = interp1( freqs, F2, sampled_freqs, 'cubic' );
%% compute the eigenrepresentation of X
% 1st coordinate of representation
tic
X1 = zeros( ( 2 * NBases + 1 ), ( 2 * NBases + 1 )  );
for k = 1:size( multi_idx, 1 )
    c_idx = multi_idx( k );
    %     for j = -NBases:NBases
    %         if ( c_idx + j >= -NBases && c_idx + j <= NBases )
    %             X1( c_idx + j + NBases + 1, k ) = f1vec( j + NBases + 1 ) + f2vec( j + NBases + 1 ) * 1i * c_idx( 1 );
    %         end
    %     end
    s_idx = c_idx + ( -NBases:NBases )';
    h_idx = intersect( find( s_idx >= -NBases ), find( s_idx <= NBases ) );
    X1( c_idx + h_idx, k ) = ( f1vec( h_idx ) + f2vec( h_idx ) * 1i * c_idx( 1 ) );
end
toc
%% zero out small entries
% X1( abs( X1 ) < 1e-10 ) = 0;

%% compute eigenvalues of different portions
if ( compute_nulls )
    tic
    [ ~, helper1 ] = spspaces( X1, 2 );
    toc
    V1 = helper1{1};
    zeroeig1 = helper1{3};
else
    tic
    [ V1, D1 ] = eig( X1, eye( size( X1 ) ) );
    toc
    zeroeig1 = [];
    for k = 1:length( D1 )
        if ( abs( D1( k, k ) ) == 0 )
            zeroeig1 = [ zeroeig1; k ];
        end
    end
end

% %% just for plotting all the elements inside of the null space
% NSamples = 5001;
% x = linspace( -pi, pi, NSamples );
% out1 = zeros( length( zeroeig1 ), NSamples );
% for m = 1:length( zeroeig1 )
%     for n = 1:NSamples
%         out1( m, n ) = out1( m, n ) + V1( :, zeroeig1( m ) )' * exp( 1i * multi_idx( : ) * x( n ) );
%     end
% end
%
% figure; hold on;
% for m = 1:length( zeroeig1 )
%     plot( x, real( squeeze( out1( m, : ) ) ) );
% end

%% generate the LP to find the distribution using Parseval's Theorem ( first set correspond to positivity requirements,
% last set correspond to real-valuedness requirement)
Aineq = zeros( 2 * NBases + 1, length( zeroeig1 ) );
Aeq = zeros( NBases + 1, length( zeroeig1 ) ); % last constraint ensures its a probability distribution
bineq = zeros( 2 * NBases + 1, 1 );
beq = zeros( NBases + 1, 1 );
beq( end ) = 1/( 2 * pi ); % last constraint ensures its a probability distribution
obj = -ones( length( zeroeig1 ), 1 ); % just check feasibility

% grab the zero index
z_idx = find( multi_idx == 0 );

% QC
lineq = zeros( length( zeroeig1 ), 1 );
lineq( : ) = -1/2 * V1( z_idx, zeroeig1 );
Q = diag( sum( V1( :, zeroeig1 ) .* V1( :, zeroeig1 ) ) );

% take care of the constant term
Aineq( 1, : ) = V1( z_idx, zeroeig1 );

% ensure positivity
r_counter = 2;
for k = 1:NBases
    
    % the 1 + cosine terms
    Aineq( r_counter, : ) = V1( z_idx, zeroeig1 );
    h_idx = multi_idx == k;
    Aineq( r_counter, : ) = Aineq( r_counter, : ) + 0.5 * V1( h_idx, zeroeig1 );
    h_idx = multi_idx == -k;
    Aineq( r_counter, : ) = Aineq( r_counter, : ) + 0.5 * V1( h_idx, zeroeig1 );
    
    r_counter = r_counter + 1;
    % the 1 + sine terms
    Aineq( r_counter, : ) = Aineq( r_counter, : ) + V1( z_idx, zeroeig1 );
    h_idx = multi_idx == k;
    Aineq( r_counter, : ) = Aineq( r_counter, : ) - 0.5 * 1i * V1( h_idx, zeroeig1 );
    h_idx = multi_idx == -k;
    Aineq( r_counter, : ) = Aineq( r_counter, : ) +  0.5 * 1i * V1( h_idx, zeroeig1 );
    
    r_counter = r_counter + 1;
end

% ensure real-valuedness
r_counter = 1;
for k = 1:NBases
    h_idx = multi_idx == k;
    Aeq( r_counter, : ) = Aeq( r_counter, : ) + V1( h_idx, zeroeig1 );
    h_idx = multi_idx == -k;
    Aeq( r_counter, : ) = Aeq( r_counter, : ) - V1( h_idx, zeroeig1 );
    
    r_counter = r_counter + 1;
end

% ensure its a probability density
Aeq( r_counter, : ) = V1( z_idx, zeroeig1 ) * ( 2 * pi );

% [ result, value, exitflag, output]  = cplexlp( obj, -Aineq, -bineq, Aeq, beq );
[ result, value, exitflag, output]  = cplexqcp( [], obj, -Aineq, -bineq, Aeq, beq, lineq, Q, 0.5 );



%% plot the element in the nullspace that comes out of the optimization
if ( exitflag == -2 )
    disp(' No feasible solution' );
elseif( exitflag == -8 )
    disp( 'Unbounded solution' );
elseif( exitflag == 1 )
    NSamples = 5001;
    x = linspace( -pi, pi, NSamples );
    out1 = zeros( 1, NSamples );
    scale = repmat( result', [ size( multi_idx, 1 ), 1 ] );
    for n = 1:NSamples
        out1( n ) = out1( n ) + sum( ( scale .* V1( :, zeroeig1 ) )' * exp( 1i * multi_idx( : ) * x( n ) ) );
    end
    figure; hold on;
    plot( x, out1( : ) );
end

