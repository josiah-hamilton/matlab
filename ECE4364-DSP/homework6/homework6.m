% Josiah Hamilton
% Homework 6
% Design a system that converts an audio signal from a sampling rate of 44.1 kHz to 48 kHz. 
% Test your system and describe its performance. 
% This is an individual assignment.
clear; close all;

% x[n] --> [^ L] --> [H(z**L)] --> [v M] --> y[p]
%               y1[l]         y2[m]
% L and M are determined using rational approximation, as suggested by matlab


fileID_Lp = fopen('h1.fcf','r');
formatSpec = '%f';
% This filter is designed for this problem, centered around the rational approximation 1/147 shown below.
% However it was expanded to reduce the order of the filter and could be more generalized.
lp = fscanf(fileID_Lp,formatSpec);

N=128;
Nfft=1024;
fs = 44100; % Sample rate
fsnew = 48000; % Desired Sample rate
f_i = 0:1/Nfft:(Nfft/2-1)/Nfft; % Normalized frequencies

x = sin(pi*0.5*(1:N)) + sin(pi*0.1*(1:N)) + sin(pi*0.33*(1:N));
X = fft(x,Nfft);

% Mathworks suggests using rational approximation

[M,L] = rat(fsnew/fs);

y1(1:M:M*length(x)) = x;
Y1 = fft(y1,Nfft);

% Design a filter such that the only non-image is preserved before downsampling.
y2 = conv(y1,lp);
Y2 = fft(y2,Nfft);

y = y1(1:M:length(y2));
Y = fft(y,Nfft);
%y_alt = sin(0.5*M/L*pi/fsnew*(1:N));
figure
subplot(431); plot(x); title("x");
subplot(432); plot(f_i,abs(X(1:Nfft/2))); 
        title("|X_u(e^j^\omega)|");
subplot(434); plot(y1); title("y1");
subplot(435); plot(f_i,abs(Y1(1:Nfft/2))); 
        title("|Y1_u(e^j^\omega)|");
subplot(437); plot(y2); title("y2");
subplot(438); plot(f_i,abs(Y2(1:Nfft/2))); 
        title("|Y2_u(e^j^\omega)|");
subplot(4,3,10); plot(y); title("y");
subplot(4,3,11); plot(f_i,abs(Y(1:Nfft/2))); 
        title("|Y_u(e^j^\omega)|");

% Test the output versus the input signal

subplot(4,3,12); plot(y-x);
        title("y-x");