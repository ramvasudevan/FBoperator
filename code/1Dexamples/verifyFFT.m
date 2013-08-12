syms x1 w1
NBases = 1000;
Fs = 4096; % increase this number to increase the condition number of the matrix (i.e. increase the accuracy of our computation)
t = -pi:( ( 2 * pi )/Fs ):( pi - ( 2 * pi)/Fs );
% t = (-pi+1e-10 ):( ( ( 2 * pi - 2e-10 ) )/Fs ):( pi - ( 2 * pi)/Fs );
% t = -3.13:6.26/Fs:( 3.13 - 6.26/Fs );
nfft = length( t );
freqs = ( -nfft/2:( nfft/2 - 1 ) ) * Fs/nfft;
sampled_freqs = -NBases:NBases;
F = matlabFunction( simplify( fourier( cos(x1), x1, w1 ) ), 'vars', w1 );
H = @(x1) cos(x1);
woot = fftshift( fft( fftshift( H( t ) ), nfft ) )/nfft * ( 2 * pi ); % CHECK THE BOOKMARK ON MY LAPTOP! NEED TO SHIFT THE SIGNAL BEFORE FFT-ing since the signal is not centered about zero which results in a shift!!!!
foo = interp1( freqs, woot, sampled_freqs, 'spline' );
boo = F( sampled_freqs + eps );
% blah = ifft( fftshift( woot ) );
