%% compute the operator 
NBases = 400;
compute_all_eigs = 1;
Xp = zeros( 2 * NBases + 1 );

for k = -NBases:NBases
    if ( k == -NBases )
        Xp( 2, k + NBases + 1 ) = 0.5 * ( k + 1 );
    elseif ( k == NBases )
        Xp( k + NBases, k + NBases + 1 ) = 0.5 * ( 1 - k );
    else
        Xp( k + NBases + 2, k + NBases + 1 ) = 0.5 * ( k + 1 );
        Xp( k + NBases, k + NBases + 1 ) = 0.5 * ( 1 - k );
    end
end

%% compute the windowed function and its DFT
Fs = 1000;
t = 0:1/Fs:( 2 * pi ); 
x = zeros( size( t ) );
x( find( t < pi/10 ) ) = 1;
nfft = 2048; 
% Take fft, padding with zeros so that length(X) is equal to nfft 
X = fft( x, nfft );
X = fftshift( X );
% Frequency vector
f = (-nfft/2:nfft/2-1)*Fs/nfft;
% Generate the plot, title and labels. 
h = figure; plot( f, real( X ) ); 
xlabel('Frequency (Hz)');

%% compute the eigenvector representation of x
xvec = zeros( 2 * NBases + 1, 1 );
fxvec = [];
for k = -NBases:NBases
    search_freq = k/( 2 * pi );
    if ( find( f == search_freq ) )
        xvec( k + NBases + 1 ) = X( find( f == search_freq ) );
    else
        min_idx = max( find( f < search_freq ) );
        max_idx = min( find( f > search_freq ) );
        min_freq = f( min_idx );
        max_freq = f( max_idx );
        t_eval = ( search_freq - min_freq )/( max_freq - min_freq );
        xvec( k + NBases + 1 ) = X( max_idx ) * t_eval + X( min_idx ) * ( 1 - t_eval );
    end
    fxvec = [ fxvec; search_freq ];
end

%% 
tau = linspace( 0, 2 * pi, length( xvec ) );
T = linspace( 0, 1, 100 );
for l = 1:length( T )
    expmXp = expm( Xp * T( l ) );
    helper = expmXp * xvec;
    out = ifft( ifftshift( helper ), 'symmetric' )/( 2 * pi );
%     out = zeros( 1, NSamples );
%     for j = 1:NSamples
%         for k = -NBases:NBases
%             out( j ) = helper( k + NBases + 1 ) * exp( 1i * k * tau( j ) );
%         end
%     end
    plot( tau, out );
    anow = axis;
    axis( [ -0.1 7 anow( 3 ) anow( 4 ) ] ) 
    filename = sprintf( 'mesa_images/image_%2.3f.pdf', T( l ) );
    print( h, '-dpdf', filename );
end




