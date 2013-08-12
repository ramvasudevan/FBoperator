%% parameters
NBases = 21;
compute_nulls = 1;

%% generate multi-idx to make clear what each element in row is
multi_idx = (-NBases:NBases)';

%% compute f_{1-4} and its DFT
tic
syms x;
% f1 = matlabFunction( fourier( heaviside( x + pi/2 ) - heaviside( x - pi/2 ) ) );
% f2 = matlabFunction( fourier( dirac( x + pi/2 ) * x ) );
% f3 = matlabFunction( fourier( -dirac( x - pi/2 ) * x ) );
% f4 = matlabFunction( fourier( ( heaviside( x + pi/2 ) - heaviside( x - pi/2 ) ) * x ) );
% f1 = matlabFunction( fourier( heaviside( x + pi/2 ) - heaviside( x - pi/2 ) ) );
% % frequency vector for analysis
% sampled_freqs = (-NBases:NBases);
% f1vec = f1( sampled_freqs + eps );
% f2vec = f2( sampled_freqs + eps );
% f3vec = f3( sampled_freqs + eps );
% f4vec = f4( sampled_freqs + eps );
f2 = matlabFunction( simplify( fourier( dirac( x + 1 ) ) ) );
f3 = matlabFunction( simplify( fourier( -dirac( x - 1 ) ) ) );
f4 = matlabFunction( simplify( fourier( ( heaviside( x + 1 ) - heaviside( x - 1 ) ) ) ) );
% f2 = @(w)(exp(1i*w));
% f3 = @(w)(-exp(-1i*w));
% f4 = @(w)(1+exp(1i*w));
% frequency vector for analysis
sampled_freqs = (-NBases:NBases);
f2vec = f2( sampled_freqs/(2 * pi) + eps );
f3vec = f3( sampled_freqs/(2 * pi) + eps );
f4vec = f4( sampled_freqs/(2 * pi) + eps );
f1vec = zeros( size( f2vec ) );
toc
%% compute the eigenrepresentation of X
% 1st coordinate of representation
tic
X= zeros( ( 2 * NBases + 1 ), ( 2 * NBases + 1 )  );
for k = 1:size( multi_idx, 1 )
    c_idx = multi_idx( k );
    %     for j = -NBases:NBases
    %         if ( c_idx + j >= -NBases && c_idx + j <= NBases )
    %             X1( c_idx + j + NBases + 1, k ) = f1vec( j + NBases + 1 ) + f2vec( j + NBases + 1 ) * 1i * c_idx( 1 );
    %         end
    %     end
    s_idx = c_idx + ( -NBases:NBases )';
    h_idx = intersect( find( s_idx >= -NBases ), find( s_idx <= NBases ) );
    X( c_idx + h_idx, k ) = ( f1vec( h_idx )+ f2vec( h_idx ) + f3vec( h_idx ) + ...
        f4vec( h_idx ) * 1i * c_idx( 1 ) );
end
toc
%% zero out small entries
X( abs( X ) < 1e-10 ) = 0;

%% compute eigenvalues of different portions
if ( compute_nulls )
    tic
    [ ~, helper1 ] = spspaces( X, 2 );
    toc
    V = helper1{1};
    zeroeig = helper1{3};
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

%% just for plotting all the elements inside of the null space
NSamples = 5001;
x = linspace( -pi, pi, NSamples );
out = zeros( length( zeroeig ), NSamples );
for m = 1:length( zeroeig )
    for n = 1:NSamples
        out( m, n ) = out( m, n ) + V( :, zeroeig( m ) )' * exp( 1i * multi_idx( : ) * x( n ) );
    end
    out( m, : ) = out( m, : )/V( NBases + 1, zeroeig( m ) );
end

figure; hold on;
for m = 1:length( zeroeig )
    plot( x, real( squeeze( out( m, : ) ) ) );
end

