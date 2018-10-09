% Homework 5 #4
% Josiah Hamilton
% adapted from Ljiljana Milic(37) Program demo2_3 

% Remarks
% Sinusoids that are higher in frequency than f_max/M are merely thrown out in the down-sampling processes
% The decimation also serves to smoothen the noise present in the down-sampled signal

clear; close all; 

% Spectrum of decimated signal
N=512;
M=5;
% Input signal 'x'
% x[n] =sin[2πf1n] + 0.9sin[2πf2n] + 0.7sin[2πf3n] + 0.8s[n]
x1 = sin(pi*0.15*(1:N)); x2 = 0.9*sin(pi*0.35*(1:N)); x3 = 0.7*sin(pi*0.55*(1:N)); noise=randn([1 512]); % Generating the signal components 
x = x1 + x2 + x3 + noise; % Original signal
X = fft(x,1024); % Computing the spectrum of the original signal
f = 0:1/512:(512-1)/512; % Normalized frequencies
figure (1)
subplot(311), plot(f,abs(X(1:512)))
ylabel('|X(e^j^\omega)|'), text(0.9,max(abs(X))/2,'(a)')

% Down-sampled signal
y = x(1:M:256); % Down-sampling
Y = fft(y,1024); % Computing the spectrum of the down-sampled signal 
subplot(312), plot(f,abs(Y(1:512)))
ylabel('|Y(e^j^\omega)|'), text(0.9,max(abs(Y))/2,'(b)')

% Decimated signal
yd = decimate(x,M); % Decimated signal
Yd = fft(yd,1024); % Computing the spectrum of the decimated signal 
subplot(313), plot(f,abs(Yd(1:512))) 
xlabel('\omega/\pi'),ylabel('|Y_d(e^j^\omega)|'), text(0.9,max(abs(Yd))/2,'(c)')