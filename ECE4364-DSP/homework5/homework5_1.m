% Homework 5 #1
% Josiah Hamilton

clear; close all; 

% general variables
N=51; 
M=4;
L=4;

% (i) sinusoidal sequence of normalized frequency 0.15
x_i = sin(2*pi*0.15*(1:N));
% (ii) sum of two sinusoidal sequences of normalized frequencies 0.1 and 0.3
x_ii = sin(2*pi*0.1*(1:N)) + sin(2*pi*0.3*(1:N));
% (iii) product of the sinusoidal sequence of normalized frequency 0.15 and the real exponential sequence {0.8n}
x_iii = sin(2*pi*0.15*(1:N)) .* 0.8.^(1:N);

% (a) Perform the factor-of-4 down-sampling. Plot the original and down-sampled sequences.
y_ia = x_i(1:M:N);
y_iia = x_ii(1:M:N);
y_iiia = x_iii(1:M:N);
N_a = length(y_ia);
% (b) Repeat part (a) for the factor-of-5 down-sampler.
M=5;
y_ib = x_i(1:M:N);
y_iib = x_ii(1:M:N);
y_iiib = x_iii(1:M:N);
N_b = length(y_ib);

figure
subplot(331)
stem(x_i); title('x_{i}')
axis([0 N -inf inf])
subplot(334)
stem(x_ii); title('x_{ii}')
axis([0 N -inf inf])
subplot(337)
stem(x_iii); title('x_{iii}')
axis([0 N -inf inf])

subplot(332)
stem(y_ia); title('y_{i•a}');
axis([0 N_a -inf inf])
subplot(335)
stem(y_iia); title('y_{ii•a}');
axis([0 N_a -inf inf])
subplot(338)
stem(y_iiia); title('y_{iii•a}');
axis([0 N_a -inf inf])

subplot(333)
stem(y_ib); title('y_{i•b}');
axis([0 N_b 0 inf])
subplot(336)
stem(y_iib); title('y_{ii•b}');
axis([0 N_b 0 inf])
subplot(339)
stem(y_iiib); title('y_{iii•b}');
axis([0 N_b 0 inf])