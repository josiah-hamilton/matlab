clear all; clc; close all; warning('off','MATLAB:plot:IgnoreImaginaryXYPart')

f = 1;
T = 3*10^8 /f;
B = 2*pi; %/ T;
z = -1:1/100:0;
Vp = 1; % Vp is given for this problem

Zo = 1; %Z is unknown, so Current will be normalized

figure('Name','Hwk4')
cases = [5,2,1/2,1/5];
Vsig=zeros(size(cases,2),size(z,2));
Gz = zeros(size(cases,2),size(z,2));
Vz = zeros(size(cases,2),size(z,2));
Iz = zeros(size(cases,2),size(z,2));
for ii=1:size(cases,2)
        
    GL(ii) = complex(GammaL(cases(ii),1));      %% For HW 4 ZL is a function of Zo
    Vm(ii,:) = VMinus(Vsig(ii,:),GL(ii));
    
    % because impedences are real for this homework
    p(ii)  = GL(ii);
    
    S(ii)  = (1+p)/(1-p);
    V(ii,:) = sqrt((Vp + GL(ii) * cos(2*B*z)).^2 + (GL(ii)*sin(2*B*z)).^2);
    I(ii,:) = Vp*sqrt((1-GL(ii) * cos(2*B*z)).^2 + (-GL(ii)*sin(2*B*z)).^2);
     
    subplot(4,1,ii);
    plot(z,V(ii,:));
    hold on;
    plot(z,I(ii,:));
    hold off;
    legend('V(z)', 'Vo*|I(z)|');
    title(strcat('Ro=', num2str(cases(ii)), '*Zo'))
    axis([z(1),z(end),0,Vp*3]);
    ylabel('V');
    xlabel('l/lambda');
    
end


function GL = GammaL(Zo, Zl)
    GL = (Zl - Zo) / (Zl + Zo);
end

function Gs = GammaS(Zs, Zo)
    Gs = (Zs - Zo)/(Zs + Zo);
end

function Vp = VPlus(Vs, Zs, Zo)
    Vp = Vs * (Zo) / (Zs + Zo);
end

function Vm = VMinus(Vp,GL)
    Vm = Vp * GL;
end

% function Gz = Gamma(z, GL, B)
%     j = sqrt(-1);
%     Gz = zeros(1,size(z,2)); %Not strictly needed, but prealloc means speed
%     Gz(1,size(z,2)) = GL * exp(j*2*B*(1:size(z,2));
% end
% 
% function Vz = Voltz(z, Vp, B, Gz, GL)
%     j = sqrt(-1);
%     Vz = zeros(1,size(z,2)); %Not strictly needed, but prealloc means speed
%     for jj=1:size(z,2)
% %         Vz(1,jj) = Vp * exp(-j*B*jj)*(1 + Gz(jj));
%           Vz(1,jj) = sqrt(Vp * );
%     end
% end
% 
% function Iz = Ampz(z, Vp, B, GL, Zo)
%     j = sqrt(-1);
%     Iz = zeros(1,size(z,2)); %Not strictly needed, but prealloc means speed
%     for jj=1:size(z,2)
%         Iz(1,jj) = Vp * exp(-j*B*jj)*(1-Gz(jj)) / Zo;
%         %Iz(1,jj) = Vp * exp(-j*B*jj) * (1 + GL) / Zo; 
%     end
% end