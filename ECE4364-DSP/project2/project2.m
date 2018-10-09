%Fixed Point Filter Implementation 
%Will Vining     Josiah Hamilton     Cris Garcia
%Adapted from Dr. Tanja Karp    
clear; clc; close all; warning('off','MATLAB:plot:IgnoreImaginaryXYPart')

%Open Filter Coefiecnts txt files

% https://www.mathworks.com/help/fixedpoint/gs/fixed-point-arithmetic_bta0vjg-1.html

fileID_Bp = fopen('bp_fil_coeff.txt','r');
fileID_Hp = fopen('hp_fil_coeff.txt','r');
fileID_Lp = fopen('lp_fil_coeff.txt','r');

formatSpec = '%f';

bp = fscanf(fileID_Bp,formatSpec);
hp = fscanf(fileID_Hp,formatSpec);
lp = fscanf(fileID_Lp,formatSpec);

fclose(fileID_Bp); fclose(fileID_Hp); fclose(fileID_Lp);

Nh=length(lp);   % filter impulse response length
Nx = 1024;      % Length of random process
bWord = 16;     % All fractional values in 
bFrac = 14;     % impulse responses and
sWord = 8;      % in random process
sFrac = 6;
baccWord = 40;
baccFrac = 35;  % Needs Justification along lines of
                % Summing up to length(h) times fractional values
                % However an all 1's test was pretty conclusive
saccWord = 20;
saccFrac = 15;

%Select Input Process
x = -1 + (1+1)*rand(1,Nx); % Uniform random
%x = zeros(1,Nx);
%x = ones(1,Nx);
%x = [1 zeros(1,Nx-2) 1];

% gains
sgl=sfi(1,8,6);
sgb=sfi(1,8,6);
sgh=sfi(1,8,6);
bgl=sfi(1,16,14);
bgb=sfi(1,16,14);
bgh=sfi(1,16,14);

Nfft=1024;  % size of fft for frequency responses
fs=8000;    % sampling frequency

%Convert to Fixed Point 
% quantize impulse responses
for ii=1:Nh
   sh_bp_fi(ii) = sfi(bp(ii),sWord,sFrac);
   sh_hp_fi(ii) = sfi(hp(ii),sWord,sFrac);
   sh_lp_fi(ii) = sfi(lp(ii),sWord,sFrac);
   bh_bp_fi(ii) = sfi(bp(ii),bWord,bFrac);
   bh_hp_fi(ii) = sfi(hp(ii),bWord,bFrac);
   bh_lp_fi(ii) = sfi(lp(ii),bWord,bFrac);
   
end

for ii=1:Nx
    sx(ii) = sfi(x(ii),sWord,sFrac);
    bx(ii) = sfi(x(ii),bWord,bFrac);
end

% Output
sy_lp = ficonv(sx,sh_lp_fi,sWord,sFrac,saccWord,saccFrac);
by_lp = ficonv(bx,bh_lp_fi,bWord,bFrac,baccWord,baccFrac);
sy_bp = ficonv(sx,sh_bp_fi,sWord,sFrac,saccWord,saccFrac);
by_bp = ficonv(bx,bh_bp_fi,bWord,bFrac,baccWord,baccFrac);
sy_hp = ficonv(sx,sh_hp_fi,sWord,sFrac,saccWord,saccFrac);
by_hp = ficonv(bx,bh_hp_fi,bWord,bFrac,baccWord,baccFrac);

y_lp  = conv(x,lp);
y_bp  = conv(x,bp);
y_hp  = conv(x,hp);

sy = sgl*sy_lp + sgb*sy_bp + sgh*sy_hp;
by = bgl*by_lp + bgb*by_bp + bgh*by_hp;

y  = y_lp + y_bp + y_hp;

%Truncate Signals Back to Nx
y_trunk = y((length(y)-length(x))/2+1:(length(y)-(length(y)-length(x))/2))

%Filter Impluse Coefficients
%Plot Impluse Coeefiecits for floating point
% graph impulse responses
n=0:Nh-1;    % time index
figure
subplot(311)
stem(n,lp)
ylabel('h_{lp}[n]')
title('filter impulse responses floating point')
subplot(312)
stem(n,bp)
ylabel('h_{bp}[n]')
subplot(313)
stem(n,hp)
ylabel('h_{hp}[n]')
xlabel('n')


%Plot Impluse Coefficients fixed point 8 bit
% graph impulse responses
n=0:Nh-1;    % time index
figure
subplot(311)
stem(n,sh_lp_fi)
ylabel('h_{lp}[n]')
title('filter impulse responses fixed point 8 bit')
subplot(312)
stem(n,sh_bp_fi)
ylabel('h_{bp}[n]')
subplot(313)
stem(n,sh_hp_fi)
ylabel('h_{hp}[n]')
xlabel('n')

%Plot Impluse Coefficients fixed point
% graph impulse responses
n=0:Nh-1;    % time index
figure
subplot(311)
stem(n,bh_lp_fi)
ylabel('h_{lp}[n]')
title('filter impulse responses fixed point 16 bit')
subplot(312)
stem(n,bh_bp_fi)
ylabel('h_{bp}[n]')
subplot(313)
stem(n,bh_hp_fi)
ylabel('h_{hp}[n]')
xlabel('n')

%Sum of impulses for respective quantum levels
n=0:Nh-1;    % time index
figure
subplot(311)
stem(n,lp+bp+hp)
ylabel('h_{all}[n]')
title('Sum of real impulse responses')
subplot(312)
stem(n,sh_lp_fi+sh_bp_fi+sh_hp_fi)
ylabel('h_{all}[n]')
title('Sum of 8-bit impulse responses')
subplot(313)
stem(n,bh_lp_fi+bh_bp_fi+bh_hp_fi)
ylabel('h_{all}[n]')
title('Sum of 16-bit impulse responses')
xlabel('n')

%Display Low-Pass Filter Outputs
figure 
subplot(311)
plot(conv(lp,x))
ylabel('h_{lp}[n]')
xlabel('n')
title('Floating Point Low-Pass')

subplot(312)
plot(sy_lp)
ylabel('h_{lp}[n]')
xlabel('n')
title('8-bit Fixed Point Low-Pass')

subplot(313)
plot(by_lp)
ylabel('h_{lp}[n]')
xlabel('n')
title('16-bit Fixed Point Low-Pass')

%Display Band-Pass Filter Outputs
figure 
subplot(311)
plot(conv(bp,x))
ylabel('h_{bp}[n]')
xlabel('n')
title('Floating Point Band-Pass')

subplot(312)
plot(sy_bp)
ylabel('h_{bp}[n]')
xlabel('n')
title('8-bit Fixed Point Band-Pass')

subplot(313)
plot(by_bp)
ylabel('h_{bp}[n]')
xlabel('n')
title('16-bit Fixed Point Band-Pass')



%Display High-Pass Filter Outputs
figure 
subplot(311)
plot(conv(hp,x))
ylabel('h_{bp}[n]')
xlabel('n')
title('Floating Point High-Pass')

subplot(312)
plot(sy_hp)
ylabel('h_{bp}[n]')
xlabel('n')
title('8-bit Fixed Point High-Pass')

subplot(313)
plot(by_hp)
ylabel('h_{bp}[n]')
xlabel('n')
title('16-bit Fixed Point High-Pass')

%Display filter combination 
figure
subplot(411)
plot(conv(lp,x)+ conv(bp,x) + conv(hp,x))
ylabel('xfl [n]')
xlabel('n')
title('Floating Point Filter Combination')

subplot(412)
plot(sy_lp + sy_bp + sy_hp)
ylabel('xfi [n]')
xlabel('n')
axis([0 1200 -1 1])
title('Fixed Piont 8-bit Filter Combination')

subplot(413)
plot(by_lp + by_bp + by_hp)
ylabel('xfi [n]')
xlabel('n')
title('Fixed Piont 16-bit Filter Combination')

subplot(414)
plot(x)
ylabel('x[n]')
xlabel('n')
title('Original Signal')

%Display Convolution 
figure
subplot(411); plot(x); title('x');
subplot(412); plot(y); title('y');
subplot(413); plot(sy); title('8-bit y');
subplot(414); plot(by); title('16-bit y');

pad_correct = round((length(y)-length(x))/2);
sy_trunk = sy(Nh+1:length(sy)-Nh);
%sy = sy_trunk;
by_trunk = by(Nh+1:length(by)-Nh);
%by = by_trunk;
Sub_sy = sy_trunk - sx;
Sub_by = by_trunk - bx;
for i=1:Nx
   
    Sub_y(i)  = (y(i+pad_correct)) - (x(i));
end

figure
subplot(311); plot(Sub_y); title('Floating point difference'); axis([0 1200 -1 1]);
subplot(312); plot(Sub_sy); title('8-bit difference'); axis([0 1200 -1 1]);
subplot(313); plot(Sub_by);  title('16-bit difference'); axis([0 1200 -1 1]);

%Convert Fixed to double for frequency response
f_fi=(0:Nfft-1)/Nfft*fs;
%Create Frequency Resoponse of 8-bit fixed Point
[sH_lp_fi,f_fl]=freqz(double(sh_lp_fi),1,Nfft,'whole',fs);
[sH_bp_fi,f_fl]=freqz(double(sh_bp_fi),1,Nfft,'whole',fs);
[sH_hp_fi,f_fl]=freqz(double(sh_hp_fi),1,Nfft,'whole',fs);

%Create Frequency Resoponse of 16-bit fixed Point
[bH_lp_fi,f_fl]=freqz(double(bh_lp_fi),1,Nfft,'whole',fs);
[bH_bp_fi,f_fl]=freqz(double(bh_bp_fi),1,Nfft,'whole',fs);
[bH_hp_fi,f_fl]=freqz(double(bh_hp_fi),1,Nfft,'whole',fs);

%Create Frequency Resoponse of Floating Point
[H_lp_fl,f_fl]=freqz(lp,1,Nfft,'whole',fs);
[H_bp_fl,f_fl]=freqz(bp,1,Nfft,'whole',fs);
[H_hp_fl,f_fl]=freqz(hp,1,Nfft,'whole',fs);

%graph frequency responses
figure

subplot(311)
plot(f_fi,20*log10(abs(sH_lp_fi)),f_fi,20*log10(abs(sH_bp_fi)),f_fi,20*log10(abs(sH_hp_fi)))
title('8-bit Fixed Point Freq Response')
legend('LP','BP','HP')
grid
xlabel('frequency [Hz]')


subplot(312)
% graph frequency responses for 8-bit fixed point
plot(f_fi,20*log10(abs(bH_lp_fi)),f_fi,20*log10(abs(bH_bp_fi)),f_fi,20*log10(abs(bH_hp_fi)))
title(' 16 -Fixed Point Freq Response')
legend('LP','BP','HP')
grid
xlabel('frequency [Hz]')
ylabel('magnitude response [dB]')

subplot(313)
% graph frequency responses for floating point
plot(f_fl,20*log10(abs(H_lp_fl)),f_fl,20*log10(abs(H_bp_fl)),f_fl,20*log10(abs(H_hp_fl)))
title('Floating Point Freq Response')
legend('LP','BP','HP')
grid
xlabel('frequency [Hz]')


%% Initialize Variables
%lambda_max=6ds;    % max lag of autocorrelation function to be considered
Nfft=1024;          % size of fft, needs to be larger than 2*lamda_max+1!
Nblock=4*128;         % blocksize for Bartlett method, needs to be less than Nfft!

w1=((0:Nfft-1)/Nfft*2*pi);  % frequency range (row vector)
w2=((0:Nx-1)/Nx*2*pi);  % frequency range for periodogram(row vector)

%Bartlett PSD 8-bit fixed point
% xblock=zeros(Nx,ceil(Nx/Nblock));
% xblock(1:Nx) = double(by_trunk);
% Sb_sx=1/ceil(Nx/Nblock)*sum((1/Nblock)*abs(fft(xblock,Nfft)).^2,2);

if Nx>=Nfft     % in this case we need to create a wrapped around version of x
    xwrap=zeros(Nfft,ceil(Nx/Nfft));
    xwrap(1:Nx)=double(y_trunk);
    xwrap=sum(xwrap,2);
    Sp=(1/Nx)*abs(fft(xwrap,Nfft)).^2;
else
    Sp=(1/Nx)*abs(fft(double(y_trunk),Nfft)).^2;
    % % calculating DTFT(creates big matrix)
    % Sp1=1/N*abs(x.'*exp(-1i*(0:N-1).'*w1)).^2;  
end 

% calculate N points in frequency domain
SpN_y=(1/Nx)*abs(fft(double(y_trunk))).^2;

% calculate Nfft points in frequency domain
if Nx>=Nfft     % in this case we need to create a wrapped around version of x
    xwrap=zeros(Nfft,ceil(Nx/Nfft));
    xwrap(1:Nx)=double(sy_trunk);
    xwrap=sum(xwrap,2);
    Sp=(1/Nx)*abs(fft(xwrap,Nfft)).^2;
else
    Sp=(1/Nx)*abs(fft(double(sy_trunk),Nfft)).^2;
    % % calculating DTFT(creates big matrix)
    % Sp1=1/N*abs(x.'*exp(-1i*(0:N-1).'*w1)).^2;  
end 

% calculate N points in frequency domain
SpN_sy=(1/Nx)*abs(fft(double(sy_trunk))).^2;

if Nx>=Nfft     % in this case we need to create a wrapped around version of x
    xwrap=zeros(Nfft,ceil(Nx/Nfft));
    xwrap(1:Nx)=double(by_trunk);
    xwrap=sum(xwrap,2);
    Sp=(1/Nx)*abs(fft(xwrap,Nfft)).^2;
else
    Sp=(1/Nx)*abs(fft(double(by_trunk),Nfft)).^2;
    % % calculating DTFT(creates big matrix)
    % Sp1=1/N*abs(x.'*exp(-1i*(0:N-1).'*w1)).^2;  
end 

% calculate N points in frequency domain
SpN_by=(1/Nx)*abs(fft(double(by_trunk))).^2;

% 
% window_Size = 256;
% N_Overlap = 128;
% 
% Pxx_float =  pwelch(y_trunk, window_Size, N_Overlap);
% Pxx_fi_s  =  pwelch(double(sy_trunk),window_Size, N_Overlap);
% Pxx_fi_b  =  pwelch(double(by_trunk),window_Size, N_Overlap);


figure
subplot(311)
plot(w2,10*log10(SpN_y))
axis([0 2*pi,-inf,inf])
title('PSD, dB floating point')
subplot(312)
plot(w2,10*log10(SpN_sy))
axis([0 2*pi,-inf,inf])
title('PSD, dB fixed point 8-bit')
subplot(313)
plot(w2,10*log10(SpN_by))
axis([0 2*pi,-inf,inf])
title('PSD, dB fixed point 16-bit')


% figure
% plot(20*log10(abs(Pxx_float)));
% hold on 
% plot(20*log10(abs(Pxx_fi_s)));
% hold on 
% plot(20*log10(abs(Pxx_fi_b)));
% hold off
% legend('Float','8-bit','16-bit')
% grid
% title('PSD')
% xlabel('w')
% ylabel('Y')

%Convolution Function
function fiout = ficonv(fixedx,fixedh,word,frac,accword,accfrac)
    convN = length(fixedx); hN = length(fixedh);
    fiout=zeros(1,convN+2*hN);
    tempx = [ zeros(1,hN), fixedx, zeros(1,hN) ];
    F = fimath('ProductMode','SpecifyPrecision','ProductWordLength',accword,'ProductFractionLength',accfrac);
    for kk = round(hN/2):(convN+hN)
        temp_accumulator = sfi(0,accword,accfrac);
        for jj = 1:hN
            temp_accumulator = accumpos(temp_accumulator , mpy(F, tempx(kk+round(hN/2)-jj), fixedh(jj)));
        end
        fiout(kk) = sfi(double(temp_accumulator),word,frac);
    end
    disp 'conv_running';
end

