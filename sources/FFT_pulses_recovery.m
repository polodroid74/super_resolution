function [ wn ] = FFT_pulses_recovery( signal, K )
%FFT_pulses_recovery Summary of this function goes here
%   Detailed explanation goes here
    N = size(signal, 2);
    
    fft_signal = abs(fftshift(fft(signal)));
    pics = find_pics(fft_signal);
    wn = (-(find_max(fft_signal, pics, K)-1)/N*2*pi+pi).'
end

function [ pics ] = find_pics( signal )
    pics = [];
    N = size(signal, 2);
    for i=2:N-1
        if (signal(i-1) < signal(i) && signal(i) > signal(i+1))
            pics = [pics i];
        end
    end
end

function [ maxs ] = find_max( signal, pics, K)
    N = size(pics, 2);
    for i=1:N
        for j=i+1:N
            if (signal(pics(i)) < signal(pics(j)))
                temp = pics(i);
                pics(i) = pics(j);
                pics(j) = temp;
            end
        end
    end
    maxs = pics(1:K);
end