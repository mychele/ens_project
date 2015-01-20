%sins detection - watch out!

% ENS homework

% First thing: import audio file, play it, time plot

%set(0,'defaulttextinterpreter','latex')


clear all;
close all;
clc;

fileName = 'segnale_107.wav';
window = 4000; %samples
from = 36*window;
samples = [from, from + window - 1];
[y, Fs] = audioread('segnale_107.wav', samples);



%fprintf('\nFor the file "segnale_107.wav": \n');
%fprintf(' sampling frequency = %d Hz \n', Fs);
%fprintf(' sample length = %5.2f s \n\n', length(y)/Fs);

% play the audio file
%fprintf('Playing the audio file \n');

%sound(y, Fs);
%pause();


%*********************************************************************
% temporal analysis

% fprintf('TEMPORAL ANALYSIS\n\n');
%
% figure(1)
% plot(y);
% title('Sound plot');
% xlabel('n');
%
% figure(2)
% t = linspace(1/Fs, length(y)/Fs, length(y));
% plot(t,y);
% axis([min(t) max(t) min(y) max(y)]);
% xlabel('nT (s)'); ylabel('y(nT)'); title('Signal in time')
% grid;


%*********************************************************************
% frequency analysis

%NFFT = 128*2^nextpow2(length(y));
Y = fft(y);%, NFFT);
f = Fs*linspace(0, 1, length(Y));

figure(3)
plot(f/1000, 20*log10(abs(Y)));
title('DFT - module in dB');%,'Interpreter','latex')

%title('DFT - module in dB, f ');
xlabel('f [kHz]');%, 'Interpreter','latex');
ylabel('|Y| [dB]');%, 'Interpreter','latex');

set(gcf, 'PaperUnits', 'points');
%set(gcf, 'PaperPosition', [0 0 1920 1080]);
set(gcf, 'PaperSize', [1200, 700]);
set(gcf, 'Color', 'w');
%format_ticks(gca);
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',33)
%export_fig instThr_c2_200.png -q101 -nocrop



%*********************************************************************
% delete sins
% the following steps are inspired by http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=910632

% as first subtract the mean of the signal from the signal
%m = mean(y);
%y = y - m;
%NFFT = 2^nextpow2(length(y));

% compute the DFT
%Y = fft(y);%, NFFT);

%some parameters
K = length(Y)/2;
p0 = 0.999; % 1 - p0 is the false detection probability
thres = 200; %sqrt(-log(1 - p0^(1/K))/log(2));

R = max(abs(Y(1:length(Y)/2))) / median(abs(Y(1:length(Y)/2)));

N_tot = ones(length(Y), 1);

notch_band = 100;

while(R > thres)
    [A, index] = max(Y(1:length(Y)/2));
    %compute the related frequency
    fm = f(index);
    disp('iteration');
    disp(fm);
    
    b1 = -2*cos(2*pi*fm/Fs);
    b = [1 b1 1];
    delta = pi*notch_band/Fs; r=1-delta;
    a2=r^2; a1 = r*b1;
    a= [1 a1 a2];
    
    
    %compute notch filter
    [N1, f1]=freqz(b, a, length(Y), 'whole', Fs);
    
    figure
    plot(f1/1000, 20*log10(abs(N1)));
    title('Notch filter - module in dB');
    xlabel('f [kHz]');
    ylabel('|H| [dB]');
    
    set(gcf, 'PaperUnits', 'points');
    %set(gcf, 'PaperPosition', [0 0 1920 1080]);
    set(gcf, 'PaperSize', [1200, 700]);
    set(gcf, 'Color', 'w');
    %format_ticks(gca);
    fig=gcf;
    set(findall(fig,'-property','FontSize'),'FontSize',33)
    %export_fig instThr_c2_200.png -q101 -nocrop
    
    pause();
    
    figure
    plot(f1/1000, angle(N1));
    title('Notch filter - phase');
    xlabel('f [kHz]');
    ylabel('\phi (H)');
    
    set(gcf, 'PaperUnits', 'points');
    %set(gcf, 'PaperPosition', [0 0 1920 1080]);
    set(gcf, 'PaperSize', [1200, 700]);
    set(gcf, 'Color', 'w');
    %format_ticks(gca);
    fig=gcf;
    set(findall(fig,'-property','FontSize'),'FontSize',33)
    %export_fig instThr_c2_200.png -q101 -nocrop
    
    pause();
    
    %filter Y
    Y = N1 .* Y;
    N_tot = N1 .* N_tot;
    
    figure
    plot(f1/1000, (abs(Y)));
    title('DFT, module (dB)');
    xlabel('f [kHz]');
    ylabel('|Y| [dB]');
    
    %compute R
    R = max(abs(Y(1:length(Y)/2))) / median(abs(Y(1:length(Y)/2)));
    % pause;
end

figure
plot(f1/1000, 20*log10(abs(N_tot)));
title('Notch filter - module in dB');
xlabel('f [kHz]');
ylabel('|H| [dB]');

set(gcf, 'PaperUnits', 'points');
%set(gcf, 'PaperPosition', [0 0 1920 1080]);
set(gcf, 'PaperSize', [1200, 700]);
set(gcf, 'Color', 'w');
%format_ticks(gca);
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',33)
%export_fig instThr_c2_200.png -q101 -nocrop

pause();

figure
plot(f1/1000, angle(N_tot));
title('Notch filter - phase');
xlabel('f [kHz]');
ylabel('\phi (H)');

set(gcf, 'PaperUnits', 'points');
%set(gcf, 'PaperPosition', [0 0 1920 1080]);
set(gcf, 'PaperSize', [1200, 700]);
set(gcf, 'Color', 'w');
%format_ticks(gca);
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',33)
%export_fig instThr_c2_200.png -q101 -nocrop

pause();


figure
plot(f/1000, 20*log10(abs(Y)));
title('DFT of filtered signal - module in dB');
xlabel('f [kHz]');
ylabel('|H| [dB]');

set(gcf, 'PaperUnits', 'points');
%set(gcf, 'PaperPosition', [0 0 1920 1080]);
set(gcf, 'PaperSize', [1200, 700]);
set(gcf, 'Color', 'w');
%format_ticks(gca);
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',33)
%export_fig instThr_c2_200.png -q101 -nocrop

pause();



out = m + ifft(Y);
%sound(out, Fs);


