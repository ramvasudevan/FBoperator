%% parameters
NBases = 51;
compute_nulls = 1;

%% generate multi-idx to make clear what each element in row is
multi_idx = (-NBases:NBases)';

%% compute f_{1-4} and its DFT
tic
syms x w;
% h = -(x^3+x^5)*(heaviside( x + 1 ) - heaviside( x - 1 ));
h = -x*exp((-x^2)*5);
d = diff( h, x );
f1 = matlabFunction( simplify( fourier( h ) ) );
f2 = matlabFunction( simplify( fourier( d ) ) );
% % frequency vector for analysis
sampled_freqs = (-NBases:NBases);
f1vec = f1( sampled_freqs );
f2vec = f2( sampled_freqs );
if ( sum( isnan( f1vec ) ) )
    foo = find( isnan( f1vec ) );
    h1 = simplify( fourier( h ) );
    for j = 1:length( foo )
        f1vec( foo( j ) ) = limit( h1, sampled_freqs( foo( j ) ) ); 
    end
end
if ( sum( isnan( f2vec ) ) )
    foo = find( isnan( f2vec ) );
    h2 = simplify( fourier( d ) );
    for j = 1:length( foo )
        f2vec( foo( j ) ) = limit( h2, sampled_freqs( foo( j ) ) ); 
    end
end
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
%     X( c_idx + h_idx, k ) = ( f1vec( h_idx )+ f2vec( h_idx ) + f3vec( h_idx ) + ...
%         f4vec( h_idx ) * 1i * c_idx( 1 ) );
    X( c_idx + h_idx, k ) = ( 1i * c_idx( 1 ) * f1vec( h_idx ) + f2vec( h_idx ));
end
toc
%% zero out small entries
% X( abs( X ) < 1e-10 ) = 0;

%% compute eigenvalues of different portions
if ( compute_nulls )
    tic
    [ ~, helper1 ] = spspaces( X, 2, 1e-5 );
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
    out( m, : ) = (2*pi)^-1 * out( m, : );
end

figure; hold on;
for m = 1:length( zeroeig )
    plot( x, real( squeeze( out( m, : ) ) ) );
end

