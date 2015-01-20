% ENS homework

% First thing: import audio file, play it, time plot

clear all;
close all;
clc;

[y, Fs] = audioread('segnale_107.wav');

%*********************************************************************
% frequency analysis

[Y, f] = pmusic(y, 6);
%f = Fs*linspace(0, 1, length(Y));

figure(3)
plot(f*Fs/2/pi/1000, Y);
title('DFT, module (dB)');
xlabel('f [kHz]');
ylabel('|Y| [dB]');

figure(4)
plot(f/1000, angle(Y));
title('DFT, phase');
xlabel('f [kHz]');
ylabel('phase(Y)');


% show the frequencies of the sin noise

mean_y = mean(abs(Y));
threshold = 3;
 
[pks,locs] = findpeaks(Y, 'threshold', threshold);

figure(5)
plot(f/1000, 20*log10(abs(Y))); hold on
title('DFT, module (dB)');
xlabel('f [kHz]');
ylabel('|Y| [dB]');
plot(f(locs)/1000,20*log10(pks+0.05),'k^','markerfacecolor',[1 0 0]), hold off

%find repetead values

freqs = ones(4, 1);
r = 1;

% first_half = f(locs);
% second_half = Fs - f(locs);
% 
% for i = 1:length(first_half)
%     for t = 1:length(second_half)
%        if first_half(i) ~= second_half(t)
%            disp('this is a sin');
%            freqs(r) = first_half(i);
%            r = r + 1;
%        end
%     end
% end
% 
figure
plot(f, abs(Y)), hold on;
plot(f(locs)*Fs/2/pi/1000, abs(Y(freqs))+0.05,'k^','markerfacecolor',[1 0 0]), hold off;
 











