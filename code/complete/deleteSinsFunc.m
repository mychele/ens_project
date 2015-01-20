function [out] = deleteSinsFunc(y, Fs)
%*********************************************************************
% delete sins
% the following steps are inspired by http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=910632

y = y .* triang(length(y));
% as first subtract the mean of the signal from the signal

NFFT = 4*2^nextpow2(length(y));

% compute the DFT
Y = fft(y);%, NFFT);
Y = Y - mean(abs(Y));
f = Fs*linspace(0, 1, length(Y));

%some parameters
% K = length(Y)/2;
% p0 = 0.999; % 1 - p0 is the false detection probability
thres = 20; %sqrt(-log(1 - p0^(1/K))/log(2));

R = max(abs(Y(1:length(Y)/2)) / median(abs(Y(1:length(Y)/2))));

[peaks, locs] = localmax(abs(Y(1:length(Y)/2)));

sinsN = 0;
while(R > thres)
    
    [maxval, index] = max(abs(peaks));
    if(index ~= 1)
        %compute the related frequency
        fm = f(locs(index));
        disp(fm);
        diffLocPlus = 20*log10(abs(maxval)) - 20*log10(abs(peaks(index + 1)))
        diffLocMinus = 20*log10(abs(maxval)) - 20*log10(abs(peaks(index - 1)))
        plot(f/1000, 20*log10(abs(Y))); %hold on;
        title('DFT, module (dB)');
        xlabel('f [kHz]');
        ylabel('|Y| [dB]');
        %plot(f(locs),20*log10(peaks+0.05),'k^','markerfacecolor',[1 0 0]), hold off
        if((21 < diffLocPlus) && (diffLocPlus < 32) && (21 < diffLocMinus) && (diffLocMinus < 32))
            sinsN = sinsN + 1;
            disp(sinsN);
            
            figure
            plot(f/1000, 20*log10(abs(Y))); hold on;
            title('DFT, module (dB)');
            xlabel('f [kHz]');
            ylabel('|Y| [dB]');
            plot(f(locs),20*log10(peaks+0.05),'k^','markerfacecolor',[1 0 0]), hold off
            
            pause();
            %compute notch filter parameters
            beta1 = cos(2*pi/Fs*fm);
            alpha1 = 0.9;
            K1 = (1 + alpha1)/2;
            a1 = [K1, -2*K1*beta1, K1];
            b1 = [1, -beta1*(alpha1 + 1), alpha1];
            
            %compute notch filter
            [N1, f1]=freqz(a1,b1,length(Y),'whole',Fs);
            figure
            plot(f1/1000, 20*log10(abs(N1)));
            title('Nothc1, module (dB)');
            xlabel('f [kHz]');
            ylabel('|N1| [dB]');
            
            %filter Y
            Y = N1 .* Y;
            [peaks, locs] = localmax(abs(Y(1:length(Y)/2)));
            
        else
            break
        end
    else
        break
    end
    %compute R
    R = max(abs(Y(1:length(Y)/2)) / median(abs(Y(1:length(Y)/2))));
    pause();
end

out = ifft(Y);