%% Name: Mrinmoy Sarkar
% email: msarkar@aggies.ncat.edu

clear all;
close all;

%% record audio and save it to a file
% Fs = 16000;
% nBits = 16;
% duration = 5;
% info = audiodevinfo;
% recorder = audiorecorder(Fs,nBits,1);
% disp('Start speaking.')
% recordblocking(recorder, duration);
% disp('End of Recording.');
% y = getaudiodata(recorder);
% audiowrite('project.wav',y,Fs);

% indxxxx = 1;
% for iii=0:0.1:13
%% read a audio file
[m_sig,fs] = audioread('project.wav');
m_sig = (m_sig');
t1 = 0:1/fs:((length(m_sig)/fs)-(1/fs));
sound(m_sig,fs)
pause(5)

%% sample signal by 4kHz sampling frequency
Fs = 4000;
k = fs/Fs;
msample_sig = downsample(m_sig,k);
Ts = 1/Fs;
t = 0:1/Fs:((length(msample_sig)/Fs)-(1/Fs));

Lfft = length(t);
Lfft = 2^ceil(log2(Lfft));
M_fre = fftshift(fft(msample_sig,Lfft));
freqm = (-Lfft/2:Lfft/2-1)/(Lfft*Ts);

%% shift the level of the signal to make positive value only
bias = abs(min(msample_sig));
mshifted_sig = msample_sig + bias; % make the signal to +ve

%% quantize
Mu = 255;
no_of_bits = 8;
M = 4;
T0 = 1; %width of the pulse
V = max(abs(mshifted_sig));
comp_sig = compand(mshifted_sig, Mu, V,'mu/compressor');
delta = (max(comp_sig)-min(comp_sig))/(2^no_of_bits-1);
[q_sig_t,q_sig] = uniform_quantize(comp_sig,no_of_bits);

%% code the signal to symbol
[code_sig,code_bin] = encode(q_sig,no_of_bits,M);
Fs1 = Fs*no_of_bits/log2(M);
t_code = 0:1/(Fs1):(length(code_sig)/Fs1 - 1/Fs1);

%% create pulse for the symbol
pulsewidth = (T0*log2(M))/(Fs*no_of_bits);
fprintf('The pulse width is : %.2f ms. \n',pulsewidth*1000);
pulse_sig = generatePulse(code_sig,T0);
Fs2 = Fs1*T0;
t_pulse = 0:1/(Fs2):(length(pulse_sig)/Fs2 - 1/Fs2);

%% psk modulation
psk_sig = pskmod(pulse_sig,M,pi/4);

%% add noise
G = 7;
sig = 0.8;
Es = mean(abs(pulse_sig).^2)*G;
No = 2*sig^2;
EsNodB = 10*log10(Es/No);
awgnchan = comm.AWGNChannel;
awgnchan.NoiseMethod='Signal to noise ratio (Es/No)';
awgnchan.EsNo = EsNodB;
noisy_sig = awgnchan(psk_sig);

%% integrate over the pulse width
int_sig = integrate(noisy_sig,T0);

%% psk demodulation 
demod_sig = pskdemod(int_sig,M,pi/4);

%% decode signal
[decode_sig,decode_bin] = decode(demod_sig,M);

%% dequantize
deq_sig = decode_sig.*delta;

%% uncompress
uncom_sig = compand(deq_sig,Mu,V,'mu/expander');

%% remove bias
unbiased_sig = uncom_sig - bias;

%% Low pass filter design
B_m = 3400; %cutoff frequency
h = fir1(256, B_m*Ts);
recv_sig = filter(h,1,unbiased_sig);
sound(recv_sig,Fs)
audiowrite('project_recv.wav',recv_sig,Fs);

Lfft_recv = length(t);
Lfft_recv = 2^ceil(log2(Lfft_recv));
recv_fre = fftshift(fft(recv_sig,Lfft_recv));
freqm_recv = (-Lfft_recv/2:Lfft_recv/2-1)/(Lfft_recv*Ts);

[number,ratio] = symerr(code_sig,demod_sig);
fprintf('No. of symbol missed %d and symbol error rate %f.\n',number,ratio);
[numberb,ratiob] = biterr(code_bin, decode_bin);
fprintf('No. of bit missed %d and bit error rate %f.\n',numberb,ratiob);


% xxx(indxxxx) = EsNodB;
% yyy(indxxxx) = ratiob;
% indxxxx = indxxxx+1;
% end
% 
% %%
% figure(100)
% semilogy(xxx(1:end-9),yyy(1:end-9),'-*');
% grid on
% ylabel('BER','FontSize', 16);
% xlabel('$\frac{E_b}{N_0}$(dB)','interpreter','latex','FontSize', 16);
% title('Bit-error rate curve');
% ax = gca;
% ax.FontSize = 16;

%% plot
figure(1)
subplot(3,1,1)
plot(t1,m_sig)
ylabel('m(t)','FontSize', 16)
xlabel('t(sec)','FontSize', 16)
title('Message signal','FontSize', 16)
ax = gca;
ax.FontSize = 16; 

subplot(3,1,2)
plot(t,msample_sig)
ylabel('m(T)','FontSize', 16)
xlabel('T(sec)','FontSize', 16)
title('Sampled message signal','FontSize', 16)
ax = gca;
ax.FontSize = 16;

subplot(3,1,3)
plot(t,comp_sig)
ylabel('m_{compressed}(T)','FontSize', 16)
xlabel('T(sec)','FontSize', 16)
title('Compressed message signal','FontSize', 16)
ax = gca;
ax.FontSize = 16;

figure(2)
subplot(3,1,1)
plot(t,q_sig_t)
ylabel('m_{quantized}(T)','FontSize', 16)
xlabel('T(sec)','FontSize', 16)
title('Quantized message signal','FontSize', 16)
ax = gca;
ax.FontSize = 16;

subplot(3,1,2)
plot(t_code,code_sig)
ylabel('m_{code}(T)','FontSize', 16)
xlabel('T(sec)','FontSize', 16)
title('Coded message signal','FontSize', 16)
ax = gca;
ax.FontSize = 16;

subplot(3,1,3)
plot(t_pulse,pulse_sig)
ylabel('m_{pulse}(T)','FontSize', 16)
xlabel('T(sec)','FontSize', 16)
title('Pulsed message signal','FontSize', 16)
ax = gca;
ax.FontSize = 16;

figure(3)
subplot(3,1,1)
plot(t_code,demod_sig)
ylabel('m_{demod}(T)','FontSize', 16)
xlabel('T(sec)','FontSize', 16)
title('Demodulated message signal','FontSize', 16)
ax = gca;
ax.FontSize = 16;

subplot(3,1,2)
plot(t,decode_sig)
ylabel('m_{decode}(T)','FontSize', 16)
xlabel('T(sec)','FontSize', 16)
title('Decoded message signal','FontSize', 16)
ax = gca;
ax.FontSize = 16;

subplot(3,1,3)
plot(t,deq_sig)
title('dequantize signal')
ylabel('m_{dequantized}(T)','FontSize', 16)
xlabel('T(sec)','FontSize', 16)
title('Dequantized message signal','FontSize', 16)
ax = gca;
ax.FontSize = 16;

figure(4)
subplot(2,1,1)
plot(t,uncom_sig)
ylabel('m_{uncompressed}(T)','FontSize', 16)
xlabel('T(sec)','FontSize', 16)
title('Uncompressed message signal','FontSize', 16)
ax = gca;
ax.FontSize = 16;

subplot(2,1,2)
plot(t,recv_sig)
ylabel('m_{recv}(T)','FontSize', 16)
xlabel('T(sec)','FontSize', 16)
title('Received message signal (after low pass filter)','FontSize', 16)
ax = gca;
ax.FontSize = 16;

figure(5)
subplot(2,1,1)
plot(freqm,abs(M_fre))
ylabel('M(f)','FontSize', 16)
xlabel('f(Hz)','FontSize', 16)
title('Spectrum of message signal','FontSize', 16)
ax = gca;
ax.FontSize = 16;

subplot(2,1,2)
plot(freqm_recv,abs(recv_fre))
ylabel('M_{recv}(f)','FontSize', 16)
xlabel('f(Hz)','FontSize', 16)
title('Spectrum of received signal','FontSize', 16)
ax = gca;
ax.FontSize = 16;


scatterplot(psk_sig)
title('4-ary PSK modulated signal','FontSize', 16)
ax = gca;
ax.FontSize = 16;
ax.Children.MarkerSize = 70;

scatterplot(noisy_sig)
title('Received noisy signal','FontSize', 16)
ax = gca;
ax.FontSize = 16;

%% helper functions

function [yt,yq]=uniform_quantize(x,no_of_bits)
    del = (max(x)-min(x))/(2^no_of_bits-1);
    yt = zeros(1,length(x));
    yq = zeros(1,length(x));
    for i=1:length(x)
        yt(i) = del*round(x(i)/del);
        yq(i) = round(x(i)/del);
    end
end
    
function [ys,yb]=encode(x,n,M)
    k = log2(M);
    ys = zeros(1,length(x)*k);
    yb = [];
    indx = 1;
    for i=1:length(x)
        bincode = dec2bin(abs(x(i)),n);
        yb = [yb de2bi(abs(x(i)),n,'left-msb')];
        for j=1:k:n
            code = bincode(j:j+k-1);
            ys(indx) =  bin2dec(code);
            indx = indx + 1;
        end
    end
end

function [y,yb]=decode(x,M)
    k = log2(M);
    n = k*M;
    y = zeros(1,length(x)/M);
    indx = 1;
    yb = [];
    for i=1:length(y)
        bincode='';
        for j=1:M
            bincode = strcat(bincode,dec2bin(x(j-1+indx),k));
        end
        indx = indx+M;
        y(i) = bin2dec(bincode);
        yb = [yb de2bi(abs(y(i)),n,'left-msb')];
    end
end

function y=generatePulse(x,k)
    y = zeros(1,length(x)*k);
    indx = 1;
    for i=1:length(x)
        for j=1:k
            y(indx) = x(i);
            indx = indx + 1;
        end
    end
end

function y=integrate(x,k)
    y = zeros(1,length(x)/k);
    indx = 1;
    for i=1:k:length(x)-k+1
        y(indx) = sum(x(i:i+k-1));
        indx = indx+1;
    end
end