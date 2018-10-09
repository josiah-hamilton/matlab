%Project 1 - Correlogram, Blackman-Tukey
%Cris Garcia, Josiah Hamilton, Will Vining, adapted from Dr. Tanja Karp
clear; close all; warning('off','MATLAB:plot:IgnoreImaginaryXYPart')

%Methods of estimating the power spectrum of a random signal from its samples
%Note, charts appear in order of creation, complementary to having been sorted by signal type.

N = 1000;
w = ((0:N-1)/(N) * 2*pi);
w2= ((0:2*N-2)/(2*N-1) * 2*pi);
% window size=2 lambda_max+1
lambda_max=25;
lambda=-lambda_max:1:lambda_max;
%Nfft=2*lambda_max+1;
Nfft=256;
w1=(0:Nfft-1)/Nfft*2*pi;

%The following sequences are available for selection
rng(1);
% real valued noise sequences
rng(1)
X(:,1)=randn(N,1);                 % Gaussian, zero mean, variance of 1
X(:,2)=(-1).^randi([0,1],[N,1]);   % -1, 1 each of probability 0.5
X(:,3)=4*rand(N,1)-2;              % uniformly distributed [-2 2)

%real valued sinusoidal sequences
A=rand(1); w0=pi*rand(1); phi=2*pi*rand(1); phi1=2*pi*rand(1);
X(:,4)=A*cos(w0*(0:N-1)+phi).'; % cosine of random amplitude, frequency and phase shift
X(:,5)=cos(pi/16*(0:N-1)+phi).'; % cosine of random phase shift w0=pi/16
X(:,6)=cos(1.1*pi/16*(0:N-1)+phi1).'; % cosine of random phase shift w0=2*pi/16

%complex valued sinusoidal sequences
X(:,7)=exp(1i*pi/16*(0:N-1)+1i*phi).'; %complex exponential of frequency w0=pi/16
X(:,8)=exp(1i*1.3*pi/16*(0:N-1)+1i*phi1).'; %complex exponential of frequency w0=1.3*pi/16

for x = X(:,[1 3 4 8])

    figure
    plot(x)
    xlabel('N')
    ylabel('x[N]')
    title('realization of random process x[N]')

%Correlogram
    % estimage autocorrelation fct first
    r_xhat_biased=xcorr(x,x,'biased'); %ifft(x.*conj(x))
    r_xhat_unbiased=xcorr(x,x,'unbiased');
    
    r_x_new_biased = conv(x,conj(x)) / N;
    %for k = 1:N-lambda_max
    %    r_x_new_biased(1:k) = (x(k+lambda_max)).*conj(x(k))/N;
    %end
    
    
    %needs our own implementation ^ v
    Sc_biased=fft([r_xhat_biased(N:end);r_xhat_biased(1:N-1)] );
    Sc_unbiased=fft([r_xhat_unbiased(N:end);r_xhat_unbiased(1:N-1)] );
    Sc_new_biased=fft([r_x_new_biased(N:end);r_x_new_biased(1:N-1)] );
    
    W_b_corr = real(sum((Sc_biased).^2)/(2*pi));
    W_u_corr = real(sum((Sc_unbiased).^2)/(2*pi));
    W_n_corr = real(sum((Sc_new_biased).^2)/(2*pi));
    
    %imag(Sc_biased)
    %imag(Sc_unbiased)
    %signalAnalyzer(r_xhat_biased) % Unavailable in online mode.
    %signalAnalyzer(r_xhat_unbiased)
    figure
    subplot(211)
    plot(-(N-1):(N-1),r_xhat_biased,-(N-1):(N-1),r_x_new_biased)
    xlabel('\lambda')
    ylabel('r_{x}[\lambda]')
    title('correlogram autocorrelation estimate')
    legend('biased correlogram biased estimation','in-house, biased correlegram estimation','location','southoutside')
    
    subplot(212)
    plot(-(N-1):(N-1),r_xhat_unbiased)
    xlabel('\lambda')
    ylabel('r_{x}[\lambda]')
    title('unbiased correlogram autocorrelation estimate')    
   
%PSD Estimates
    figure
    subplot(211)
    plot(w2,Sc_biased,w2,Sc_new_biased)
    legend(strcat('biased ',W_u_corr),strcat('new biased ',W_n_corr),'location','southoutside')
    title('PSD linear, Correlogram method');
    xlabel('\omega')
    ylabel('S_x(e^{j\omega})')
    
    subplot(212)
    plot(w2,Sc_unbiased)
    legend(strcat('unbiased ',W_u_corr),'location','southoutside')
    title('PSD linear, Correlogram method');
    xlabel('\omega')
    ylabel('S_x(e^{j\omega})')
    
    figure
    subplot(211)
    plot(w2,10*log10([Sc_biased,Sc_new_biased]))
    legend('biased','new biased','location','southoutside')
    title('PSD dB, Correlogram method');
    xlabel('\omega')
    ylabel('S_x(e^{j\omega})')
    
    subplot(212)
    plot(w2,10*log10(Sc_unbiased))
    legend('unbiased','location','southoutside')
    title('PSD dB, Correlogram method');
    xlabel('\omega')
    ylabel('S_x(e^{j\omega})')

%Several POVs for Blackman-Tukey
    for lambda_max = [10 50 200]
        lambda=-lambda_max:1:lambda_max;

%Blackman-Tukey
%Don't use windows provided in class and instead use Blackman and Tukey windows with a few parameters
%Biased
        % estimate autocorrelation fct (biased) up to value lambda_max
        [rBT_rectwin,biasLag]=xcorr(x,x,lambda_max,'biased');
        
        % create autocorrelation vector of length Nfft with lambda=0 in first value
        rBT_Nfft=zeros(Nfft,1);
        rBT_Nfft(1:lambda_max+1)=rBT_rectwin(lambda_max+1:2*lambda_max+1);
        rBT_Nfft(Nfft-lambda_max+1:Nfft)=rBT_rectwin(1:lambda_max);
        %PSD calculation via fft
        SBT_rectwin=fft(rBT_Nfft,Nfft);
        
        % blackman window
        %wvtool(blackman(2*lambda_max+1));
        rBT_blackman = rBT_rectwin.*blackman(2*lambda_max+1);
        % create autocorrelation vector of length Nfft with lambda=0 in first value
        rBT_Nfft=zeros(Nfft,1);
        rBT_Nfft(1:lambda_max+1)=rBT_blackman(lambda_max+1:2*lambda_max+1);
        rBT_Nfft(Nfft-lambda_max+1:Nfft)=rBT_blackman(1:lambda_max);
        %PSD calculation via fft
        SBT_blackman=fft(rBT_Nfft,Nfft);
        
        % tukey window 0.2
        %wvtool(tukeywin(2*lambda_max+1,0.2));
        rBT_tukey2 = rBT_rectwin.*tukeywin(2*lambda_max+1,0.2);
        % create autocorrelation vector of length Nfft with lambda=0 in first value
        rBT_Nfft=zeros(Nfft,1);
        rBT_Nfft(1:lambda_max+1)=rBT_tukey2(lambda_max+1:2*lambda_max+1);
        rBT_Nfft(Nfft-lambda_max+1:Nfft)=rBT_tukey2(1:lambda_max);
        %PSD calculation via fft
        SBT_tukey2=fft(rBT_Nfft,Nfft);
        
        % tukey window 0.5
        %wvtool(tukeywin(2*lambda_max+1,0.5));
        rBT_tukey5 = rBT_rectwin.*tukeywin(2*lambda_max+1,0.5);
        % create autocorrelation vector of length Nfft with lambda=0 in first value
        rBT_Nfft=zeros(Nfft,1);
        rBT_Nfft(1:lambda_max+1)=rBT_tukey5(lambda_max+1:2*lambda_max+1);
        rBT_Nfft(Nfft-lambda_max+1:Nfft)=rBT_tukey5(1:lambda_max);
        %PSD calculation via fft
        SBT_tukey5=fft(rBT_Nfft,Nfft);
        
        % tukey window 0.8
        %wvtool(tukeywin(2*lambda_max+1,0.8));
        rBT_tukey8 = rBT_rectwin.*tukeywin(2*lambda_max+1,0.8);
        % create autocorrelation vector of length Nfft with lambda=0 in first value
        rBT_Nfft=zeros(Nfft,1);
        rBT_Nfft(1:lambda_max+1)=rBT_tukey8(lambda_max+1:2*lambda_max+1);
        rBT_Nfft(Nfft-lambda_max+1:Nfft)=rBT_tukey8(1:lambda_max);
        %PSD calculation via fft
        SBT_tukey8=fft(rBT_Nfft,Nfft);
%Power Estimates
        W_bl_bt = real(sum((SBT_blackman).^2)/(2*pi));
        W_t2_bt = real(sum((SBT_tukey2).^2)/(2*pi));
        W_t5_bt = real(sum((SBT_tukey5).^2)/(2*pi));
        W_t8_bt = real(sum((SBT_tukey8).^2)/(2*pi));
%Autocorrelation Functions
        figure
        subplot(211)
        plot(lambda,rBT_blackman, lambda, rBT_tukey2, lambda, rBT_tukey5, lambda, rBT_tukey8)
        legend('blackman', 'tukey2', 'tukey5', 'tukey8','location','southoutside')
        xlabel('\lambda')
        ylabel('r_{x}[\lambda]')
        title('Biased blackman-tukey autocorrelation estimate')
%Graph PSD Estimates
        figure
        subplot(211)
        plot(w1,[SBT_blackman,SBT_tukey2,SBT_tukey5,SBT_tukey8])
        axis([0 2*pi,-inf,inf])
        legend(strcat('BT-blackman ',W_bl_bt),strcat('BT-tukey2 ',W_t2_bt),strcat('BT-tukey5 ',W_t5_bt),strcat('BT-tukey8 ',W_t8_bt),'location','southoutside')
        title('PSD linear, Biased Blackman-Tukey method');
        xlabel('\omega')
        ylabel('S_x(e^{j\omega})')
        
        subplot(212)
        plot(w1,10*log10([SBT_blackman,SBT_tukey2,SBT_tukey5,SBT_tukey8]))
        axis([0 2*pi,-inf,inf])
        legend('BT-blackman','BT-tukey2','BT-tukey5','BT-tukey8','location','southoutside')
        title('PSD dB. Biased Blackman-Tukey method')
        xlabel('\omega')
        ylabel('S_x(e^{j\omega})')

    end
end

