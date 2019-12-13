%% Name: Mrinmoy Sarkar
% email: msarkar@aggies.ncat.edu
clear all;
close all;


%% record audio and save it to a file
% Fs = 8000;
% nBits = 8;
% duration = 5;
% info = audiodevinfo;
% recorder = audiorecorder(Fs,nBits,1);
% disp('Start speaking.')
% recordblocking(recorder, duration);
% disp('End of Recording.');
% y = getaudiodata(recorder);
% audiowrite('assign1.wav',y,Fs);

%% read a audio file
[m_sig,Fs] = audioread('assign1.wav');
m_sig = m_sig';
sound(m_sig,Fs)
pause(10)

fc = Fs/2;
T = length(m_sig)/Fs;
Ts = 1/Fs;
t = 0:Ts:T-Ts;

Lfft = length(t);
Lfft = 2^ceil(log2(Lfft));
M_fre = fftshift(fft(m_sig,Lfft));
freqm = (-Lfft/2:Lfft/2-1)/(Lfft*Ts);

%% Low pass filter design
B_m = 4000;
h = fir1(256, B_m*Ts);

%% AM modulation
mu = 0.7;
mp = max(abs(m_sig));
A = mp/mu;
s_am = (A+m_sig).*cos(2*pi*fc*t);
Lfft = length(t);
Lfft = 2^ceil(log2(Lfft)+1);
S_am = fftshift(fft(s_am,Lfft));
freq_am = (-Lfft/2:Lfft/2-1)/(Lfft*Ts);

%% AM demodulation
s_dem_am = s_am.*(s_am>0);
s_rec_am = filter(h,1,s_dem_am);
S_rec_am = fftshift(fft(s_rec_am,Lfft));

sound(s_rec_am,Fs)
pause(10)
audiowrite('assign1_rec_am.wav',s_rec_am,Fs);

%% FM modulation
beta = 0.07;
fDev = B_m*beta;
s_fm = fmmod(m_sig,fc,Fs,fDev);
Lfft = length(t); 
Lfft = 2^ceil(log2(Lfft)+1);
S_fm = fftshift(fft(s_fm,Lfft));
freq_fm = (-Lfft/2:Lfft/2-1)/(Lfft*Ts);

%% FM demodulation
s_fmdem = fmdemod(s_fm,fc,Fs,fDev);
s_rec_fm = filter(h,1,s_fmdem);
S_rec_fm = fftshift(fft(s_rec_fm,Lfft));

sound(s_rec_fm,Fs)
pause(10)
audiowrite('assign1_rec_fm.wav',s_rec_fm,Fs);

%% plot of message signal, transmitted and received am, fm signal and their spectrum
figure(1)
Frange = [-4000 4000 0 150];


subplot(5,2,1);
plot(t,m_sig);
ylabel('m(t)','FontSize', 16)
xlabel('t(sec)','FontSize', 16)
title('Message signal','FontSize', 16)

subplot(5,2,2);
plot(freqm,abs(M_fre));
axis(Frange)
ylabel('M(f)','FontSize', 16)
xlabel('f(Hz)','FontSize', 16)
title('Message spectrum','FontSize', 16)


subplot(5,2,3);
plot(t,s_am);
ylabel('s_{am}(t)','FontSize', 16)
xlabel('t(sec)','FontSize', 16)
title('AM transmitted signal','FontSize', 16)

subplot(5,2,4);
plot(freq_am,abs(S_am));
axis(Frange)
ylabel('S_{am}(f)','FontSize', 16)
xlabel('f(Hz)','FontSize', 16)
title('AM transmitted spectrum','FontSize', 16)


subplot(5,2,5);
plot(t,s_rec_am);
ylabel('s_{rec\_am}(t)','FontSize', 16)
xlabel('t(sec)','FontSize', 16)
title('AM received signal','FontSize', 16)

subplot(5,2,6);
plot(freq_am,abs(S_rec_am));
axis(Frange)
ylabel('S_{rec\_am}(f)','FontSize', 16)
xlabel('f(Hz)','FontSize', 16)
title('AM received spectrum','FontSize', 16)


subplot(5,2,7);
plot(t,s_fm);
ylabel('s_{fm}(t)','FontSize', 16)
xlabel('t(sec)','FontSize', 16)
title('FM transmitted signal','FontSize', 16)

subplot(5,2,8);
plot(freq_fm,abs(S_fm));
axis(Frange)
ylabel('S_{fm}(f)','FontSize', 16)
xlabel('f(Hz)','FontSize', 16)
title('FM transmitted spectrum','FontSize', 16)


subplot(5,2,9);
plot(t,s_rec_fm);
ylabel('s_{rec\_fm}(t)','FontSize', 16)
xlabel('t(sec)','FontSize', 16)
title('FM received signal','FontSize', 16)

subplot(5,2,10);
plot(freq_fm,abs(S_rec_fm));
axis(Frange)
ylabel('S_{rec\_fm}(f)','FontSize', 16)
xlabel('f(Hz)','FontSize', 16)
title('FM received spectrum','FontSize', 16)


