Fs = 1000; % Sampling frequency
t = 0:1/Fs:1; % Time vector of 1 second 
f = 5; % Create a sine wave of f Hz.
x = sin(2*pi*t*f); 
nfft = 2048; % Length of FFT
% Take fft, padding with zeros so that length(X) is equal to nfft 
X = fft(x,nfft); 
% FFT is symmetric, throw away second half
% X = X(1:nfft/2); 
X = fftshift(X);
% Take the magnitude of fft of x
mx = abs(X);
% Frequency vector
f = (-nfft/2:nfft/2-1)*Fs/nfft;
% f = (0:nfft/2-1)*Fs/nfft; 
% Generate the plot, title and labels. 
figure(1);
plot(t,x);
title('Sine Wave Signal'); 
xlabel('Time (s)'); 
ylabel('Amplitude'); 
figure(2);
plot(f,mx); 
title('Power Spectrum of a Sine Wave'); 
xlabel('Frequency (Hz)'); 
ylabel('Power');