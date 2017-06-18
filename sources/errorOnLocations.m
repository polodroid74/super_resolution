function [ errorWrtSNR ] = errorOnLocations( method, wn, input, K, SNRvalues )
    wn = sort(wn, 'ascend');
    N= size(input, 2);
    inputAmplitude = rms(input);
    errorWrtSNR = zeros(size(SNRvalues));
    i = 1;
    NumberOfIterations = 20;
    for SNR = SNRvalues % SNR in dB
        estimtk = zeros(K, 1);
        for iter=1:NumberOfIterations
            % Display for locations
            ecart_type = inputAmplitude/10^(SNR/20);
            noisyInput = input + ecart_type^2 * rand(1, N);
            if method==0
                res = FFT_pulses_recovery(noisyInput, K);
                % Base Prony
            elseif method==1
                res = Prony(noisyInput, K);
                % Prony TLS
            elseif method==2
                res = PronyTLS(noisyInput, K);
                % Yule-Walker
            elseif method==3
                res = YuleWalker(noisyInput, K);
                % Pisarenko
            elseif method==4
                res = Pisarenko(noisyInput, K);
                % Music
            elseif method==5
                res = -rootmusic(noisyInput, K);
                % Esprit
            elseif method==6
                res = Esprit(noisyInput, K);
                % MatrixPencil
            elseif method==7
                res = MatrixPencil(noisyInput, K, N/3);
            end
            estimtk = estimtk + sort(res, 'ascend');
        end
        estimtk = estimtk/NumberOfIterations;
        errorWrtSNR(i) = sqrt((wn - real(estimtk)')*((wn - real(estimtk)')'));
        i = i+1;
    end
    %plot(SNRvalues, errorLocationsProny);
end

