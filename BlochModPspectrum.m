

%%  Declarations

clc;
clear;
close all;
 
delta = 0;%2*pi*4e9; %Detuning in rad/sec (2pi*Hz)
t2=13e-9 ;%500e-8; %Relaxation Time T2 in s,
t1=t2; %in s
E0=1e6;   %Electric Field Strength V/m
f=4.70631*1e14; %frequency of Transition in Hz, Corrresponds to 1.9ev Transition
y0 = [0 ; -1]; %initial conditions origianlly [0,-1]
ft=linspace(-1e-8,1e-8,10000);
h = 6.62607015e-34;
hbar = h/(2*pi);
dipoletransition = 5.2 * 3.33564095e-30; % 5.2 debye transition dipole moment to si units C.m

%ft=linspace(-10e-10,10e-10,1000); 

%% Rabi Frequency Calc
rabifreq=2*(dipoletransition*E0/hbar)/(2*pi); %in HZ;
Equivalentfreq = sqrt(rabifreq^2 + ((delta)/(2*pi))^2); %in Hz

%% 
efieldnormal = E0*exp(-1j*2*pi*f*ft)+E0*exp(1j*2*pi*f*ft);
%efieldgaussian = E0*exp(-1j*2*pi*f*ft).*exp((-(ft-1e-10).^2)./2e-20);

%E0*exp(1j*2*pi*f*ft).*exp((-((ft-1e-10).^2)./2e-21));

%% E field Curves

figure;
subplot(2,1,1)
plot(ft,efieldnormal);
title('Electric Field Amplitude');
xlabel('t in s');

% subplot(2,1,2)
% plot(ft,abs(efieldgaussian));
% title('Electric Field Amplitude with Gaussuian Envelope');
% xlabel('t in s');

%% Dipole plot for both efield amplitude

opts = odeset('Refine',3000); %ODE Tolerence

 
tic
[ta,pw1]=ode89(@(ta,pw) pol(ta,pw,delta,t2,t1,E0,f),[0 30e-9], y0,opts);
toc

pw1center = (pw1(:,1)).*exp(-1i*2*pi*f*ta).*1 + (conj((pw1(:,1))).*exp(1i*2*pi*f*ta).*1);


figure;
subplot(2,1,1)
% yyaxis left

plot(ta,pw1center);
title('Dipole Magnitude');
ylabel('Dipole Magnitude in C.m');
xlabel('t in s');

subplot(2,1,2)
title('Population');
plot(ta,pw1(:,2));
ylabel('w');
xlabel('t in s');

%save('pw1center.mat','pw1center','ta',"rabifreq",'E0',"pw1",'-v7.3');

figure;
plot(ta,real(pw1(:,1)));

t = ta;
v = pw1(:,1);
f_fund=4.70631*1e14;
Fs=4e15;
targetSampleRate = Fs;
[vr,ty] = resample(v,ta,targetSampleRate,'linear');
vrcenter = ((vr).*exp(-1i*2*pi*f_fund*ty).*1) + (conj(vr).*exp(1i*2*pi*f_fund*ty).*1);
z=sin(2*pi*3e14*ty);

plot(ty,vrcenter);
title('Dipole Magnitude');
xlabel("t (s)");
ylabel("Amplitude in C.m ")
%xlim([0 0.2e-13]); 

Lvr=length(vrcenter);
Fn = Fs/2;                                                  % Nyquist Frequency
FTvr = abs(fft(vrcenter))/Lvr;                                           % Fourier Transform

Fv = linspace(0, 1, fix(Lvr/2)+1)*Fn;                         % Frequency Vector
Iv = 1:length(Fv);                                          % Index Vector
figure(2)
plot(Fv, (FTvr(Iv))*2)
grid
xlabel("f (Hz)");
ylabel("Amplitude in C.m ")

[pks,locs]=findpeaks(FTvr(1:Lvr/2+1),Fv);
[Xsorted,I] = sort(pks,'descend');
Ysorted = locs(I);

display(Ysorted(1)-Ysorted(2));


figure;
Y = fft(vrcenter);

P2 = abs(Y/Lvr);
P1 = P2(1:Lvr/2+1);
P1(2:end-1) = 2*P1(2:end-1);


f = Fs*(0:(Lvr/2))/Lvr;
plot(f,P1) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')
xlim([4.702e14 4.72e14]);

ydft = fftshift(Y);
figure;
df = Fs/Lvr;
freqvec = -Fs/2+df:df:Fs/2;
plot(freqvec,abs(ydft))
xlabel('Hz');
