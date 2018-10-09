% Homework 5 #5
% Josiah Hamilton
% adapted from Ljiljana Milic(39) Program demo2_4

% Remarks
% The two pairs of positive image sinusoids that appear in the up-sampling process do not have the same amplitudes as the original sinusoids

clear; close all; 

N=512;
L=5;

% Spectrum of interpolated signal
% Input signal ‘x’
%x[n] = sin[2π0.05n] + 0.9sin[2π0.0.9n] + 0.8s[n]
x1 = sin(2*pi*0.05*(1:N)); x2 = 0.9*sin(2*pi*0.09*(1:N)); noise=randn([1 512]); % Generating the signal components 
x = x1 + x2 + noise; % Original signal
X = fft(x,1024); % Computing the spectrum of the original signal
f_i = 0:1/1024:(512-1)/1024; % Normalized frequencies
figure (1)
subplot(311), plot(f_i,abs(X(1:512)))
ylabel('|X(e^j^\omega)|'), text(0.9,max(abs(X)),'(a)')

L = 5; % Up-sampling factor
xu = zeros(1,L*length(x));
xu([1:L:length(xu)]) = x; % Up-sampled signal
Xu = fft(xu,1024); % Computing the spectrum of the up-sampled signal 
subplot(312), plot(f_i,abs(Xu(1:512)))
ylabel('|X_u(e^j^\omega)|'), text(0.9,max(abs(Xu)),'(b)')

y = interp(x,L); % Interpolated signal
Y = fft(y,1024); % Computing the spectrum of the up-sampled signal 
subplot(313), plot(f_i,abs(Y(1:512)))
ylabel('|(Y(e^j^\omega)|'), xlabel('\omega/( \pi)'), text(0.9,max(abs(Y)),'(c)')