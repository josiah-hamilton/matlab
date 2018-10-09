%Homework 3
%Cris Garcia, Josiah Hamilton, William Vining
clear; close all;

%Initialize Variables
load("noisy_sin.mat");
load("FilterOutput.mat")

%Power - Integrate the PSD estimates and divide by (2*pi).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
subplot(4,1,1)
plot(x1)
xlabel('n')
ylabel('x1[n]')

subplot(4,1,2)
plot(x2)
xlabel('n')
ylabel('x2[n]')

subplot(4,1,3)
plot(x3)
xlabel('n')
ylabel('x3[n]')

subplot(4,1,4)
plot(x4)
xlabel('n')
ylabel('x4[n]')


N = 1000;
w=((0:N-1)/N*2*pi);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
psd1 = (1/N)*(abs(fft(x1))).^2;
psd2 = (1/N)*(abs(fft(x2))).^2;
psd3 = (1/N)*(abs(fft(x3))).^2;
psd4 = (1/N)*(abs(fft(x4))).^2;

[maxvalue1,indexmax1]=max(abs(fft(x1-mean(x1))));
[maxvalue2,indexmax2]=max(abs(fft(x2-mean(x2))));
[maxvalue3,indexmax3]=max(abs(fft(x3-mean(x3))));
[maxvalue4,indexmax4]=max(abs(fft(x4-mean(x4))));


figure
subplot(4,1,1)
plot(psd1)
title('x1')


subplot(4,1,2)
plot(psd2)
title('x2')


subplot(4,1,3)
plot(psd3)
title('x3')

subplot(4,1,4)
plot(psd4)
title('x4')



figure
subplot(4,1,1)
plot(psd1)
title('x1 db')

subplot(4,1,2)
plot(10*log10(psd2))
title('x2 db')

subplot(4,1,3)
plot(10*log10(psd3))
title('x3 db')

subplot(4,1,4)
plot(10*log10(psd4))
title('x4 db')




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
subplot(3,1,1)
plot(y1)
title('y1')

subplot(3,1,2)
plot(y2)
title('y2')

subplot(3,1,3)
plot(y3)
title('y3')


figure
subplot(3,1,1)
plot((1/N)*(abs(fft(y1))).^2)
title('y1')

subplot(3,1,2)
plot((1/N)*(abs(fft(y2))).^2)
title('y2')

subplot(3,1,3)
plot((1/N)*(abs(fft(y3))).^2)
title('y3')






