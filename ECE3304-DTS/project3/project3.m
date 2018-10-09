%Josiah Hamilton
%Ali Louangrath
%Project 3
%04/30/2016

clear all; close;
addpath(genpath(pwd));
load('num1.2.mat');
load('Num2.mat');
load('Num3.mat');

[x1,fs1] = audioread('noisy1.wav');
[x2,fs2] = audioread('noisy2.wav');
[x3,fs3] = audioread('noisy3.wav');


samples1 = length(x1);   %Use MATLAB to find the length of the sound clip

f1 = [ -(samples1 / 2) + 1:1:samples1 / 2 ] * fs1 / samples1;
x1 = x1/max(abs(x1));     %Normalize the values of x


fftx1 = abs(fftshift(fft(x1)))/samples1;
fftx1 = fftx1/max(fftx1);                           
plot(f1,fftx1);
title('Fast Fourier Transform')
xlabel('time')
ylabel('Intesity')
pause(8);

spectrogram(x1);
soundsc(x1,fs1);
pause(8);

a1 = 1;
b1 = Num; 

y1 = filtfilt(b1,a1,x1);

samples1_1 = length(y1);   %Use MATLAB to find the length of the sound clip
f1_1 = [ -(samples1_1 / 2) + 1:1:samples1_1 / 2 ] * fs1 / samples1_1;
y1 = y1/max(abs(y1));     %Normalize the values of x

ffty1 = abs(fftshift(fft(y1)))/samples1_1;
ffty1 = ffty1/max(ffty1);                           
plot(f1_1,ffty1);
title('Fast Fourier Transform');
xlabel('time');
ylabel('Intesity');
pause(8);


spectrogram(y1);
sound(y1,fs1);

pause(8);

samples2 = length(x2);   %Use MATLAB to find the length of the sound clip

f2 = [ -(samples2 / 2) + 1:1:samples2 / 2 ] * fs2 / samples2;
x2 = x2/max(abs(x2));     %Normalize the values of x

fftx2 = abs(fftshift(fft(x2)))/samples2;
fftx2 = fftx2/max(fftx2);                           
plot(f2,fftx2);
title('Fast Fourier Transform');
xlabel('time');
ylabel('Intesity');
pause(8);


spectrogram(x2);
soundsc(x2,fs2);
pause(8);

a2 = 1;
b2 = Num2; 

y2 = filtfilt(b2,a2,x2);
newy2 = y2 - y1; 

samples2_1 = length(newy2);   %Use MATLAB to find the length of the sound clip

f2_1 = [ -(samples2_1 / 2) + 1:1:samples2_1 / 2 ] * fs2 / samples2_1;
newy2 = newy2/max(abs(newy2));     %Normalize the values of x

ffty2 = abs(fftshift(fft(newy2)))/samples2_1;
ffty2 = ffty2/max(ffty2);                           
plot(f2_1,ffty2);
title('Fast Fourier Transform');
xlabel('time');
ylabel('Intesity');
pause(8);

spectrogram(newy2);
sound(newy2,fs2);

pause(8);


samples3 = length(x3);   %Use MATLAB to find the length of the sound clip

f3 = [ -(samples3 / 2) + 1:1:samples3 / 2 ] * fs3 / samples3;
x3 = x3/max(abs(x3));     %Normalize the values of x

fftx3 = abs(fftshift(fft(x3)))/samples3;
fftx3 = fftx3/max(fftx3);                           
plot(f3,fftx3);
title('Fast Fourier Transform')
xlabel('time')
ylabel('Intesity')
pause(8);

spectrogram(x3);
soundsc(x3,fs3);
pause(8);

a3 = 1;
b3 = Num3;

y3 = filtfilt(b3,a3,x3);
newy3 = y3 - y1 - y2;

samples3_1 = length(newy3);   %Use MATLAB to find the length of the sound clip

f3_1 = [ -(samples3_1 / 2) + 1:1:samples3_1 / 2 ] * fs3 / samples3_1;
newy3 = newy3/max(abs(newy3));     %Normalize the values of x

ffty3 = abs(fftshift(fft(newy3)))/samples3_1;
ffty3 = ffty3/max(ffty3);                           
plot(f3_1,ffty3);
title('Fast Fourier Transform');
xlabel('time');
ylabel('Intesity');
pause(8);

spectrogram(newy3);
sound(newy3,fs3);

pause(8);


