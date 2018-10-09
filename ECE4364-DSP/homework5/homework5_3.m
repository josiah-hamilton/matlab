% Homework 5 #3
% Josiah Hamilton

% Remarks
% Aside from the zero padding observed in y_1 through y_5, it is observed how they can sum back together to form the input signal again.

clear; close all; 

% general variables
N=51; 
M=4;
L=4;
% composed as a sum of two sinusoidal sequences of normalized frequencies 0.0625 and 0.2.
x = sin(2*pi*0.0625*(1:N)) + sin(2*pi*0.2*(1:N));

% down-sample the following sequences: {x[n]}, {x[n−1]}, {x[n−2]}, {x[n−3]}, {x[n−4]}, {x[n−5]}.
shift = 5;
xshift = [ zeros(1,shift) x ];
y_0 = xshift(shift+1:M:N);
y_1 = xshift(shift+0:M:N);
y_2 = xshift(shift-1:M:N);
y_3 = xshift(shift-2:M:N);
y_4 = xshift(shift-3:M:N);
y_5 = xshift(shift-4:M:N);

figure
subplot(331)
stem(x); title('x');
axis([0 N -inf inf])
subplot(334)
stem(y_0); title('y_{n}');
axis([0 length(y_0) -inf inf])
subplot(335)
stem(y_1); title('y_{n-1}');
axis([0 length(y_1) -inf inf])
subplot(336)
stem(y_2); title('y_{n-2}');
axis([0 length(y_2) -inf inf])
subplot(337)
stem(y_3); title('y_{n-3}');
axis([0 length(y_3) -inf inf])
subplot(338)
stem(y_4); title('y_{n-4}');
axis([0 length(y_4) -inf inf])
subplot(339)
stem(y_5); title('y_{n-5}');
axis([0 length(y_5) -inf inf])
