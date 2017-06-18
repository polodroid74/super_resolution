function main

noise=0; %1: Gaussian noise addition, 0: No additional noise
N=255; %Number of mesurements, ie line pixels
M=(N-1)/2;
gridX=linspace(-M,M,N);
wn = [0.7 1.3 2.6];
an = [1.0 1.0 1.0];
K=size(wn, 2);
input = an * exp(-1i*wn'*gridX);
%plot(gridX, input, 'r');

method=7;   % Method choice :
% 0: FFT : pour motivations
% 1: Base prony
% 2: Prony Total Least Squares (TLS)
% 3: Yule-Walker
% 4: Pisarenko
% 5: Music
% 6: Esprit
% 7: Matrix-Pencil

displayErrorGraph=1; % 1: Plot error with regards to SNR and displays the CramerRao bound, 0: No plot

inputNoise = input;
% Noise addition
if noise==1
        inputNoise=input + 0.05*randn(1, N);
end

%% Arguments reconstruction

% Pour motivation : fft
if method==0
    estimtk = FFT_pulses_recovery(inputNoise, K);
    
% Base Prony
elseif method==1
    estimtk = Prony(inputNoise, K);

% Prony TLS
elseif method==2
    estimtk = PronyTLS(inputNoise, K);

% Yule-Walker
elseif method==3
    estimtk = YuleWalker(inputNoise, K);

% Pisarenko
elseif method==4
    estimtk = Pisarenko(inputNoise, K);

% Music
elseif method==5
    %estimtk = -Music(input, K);
    estimtk = -rootmusic(inputNoise, K);

% Esprit
elseif method==6
    estimtk = Esprit(inputNoise, K);

% Matrix-Pencil
elseif method==7
    estimtk = MatrixPencil(inputNoise, K, N/3);
       
end


%% Amplitude reconstruction
Umatrix = ones(2*M + 1, K);
for l=-M:M
    for c=1:K
        Umatrix(l+M+1, c) = exp(-1i*l*estimtk(c));
    end
end
UHmatrix = Umatrix';
estimak = ((UHmatrix*Umatrix)\UHmatrix)*transpose(input);
estimak = abs(real(estimak));


%% Output : 
disp('Signal estimated components :')
for k=1:K
    fprintf('K=%d : Ampiltude = %f\n      Argument  = %f\n', k, estimak(k), estimtk(k));
end
%disp('Arguments :');
%disp(estimtk)
%disp('Amplitudes :');
%disp(estimak)


%% Graph display - Error with regards to SNR and Cramer-Rao bound display
if displayErrorGraph == 1
    clf();
    NumberOfSNRValues = 80;
    SNRvalues = linspace(1, 80, NumberOfSNRValues);
    figure(1);
    hold on;
    
    %%%%%%%%%%% Error with different methods %%%%%%%%    
    %plot(SNRvalues, log10(errorOnLocations(0, wn, input, K, SNRvalues)), 'y'); %FFT
    plot(SNRvalues, log10(errorOnLocations(1, wn, input, K, SNRvalues)), 'r'); %Prony
    plot(SNRvalues, log10(errorOnLocations(2, wn, input, K, SNRvalues)), 'g'); %PronyTLS
    plot(SNRvalues, log10(errorOnLocations(3, wn, input, K, SNRvalues)), 'k'); %YuleWalker
    plot(SNRvalues, log10(errorOnLocations(4, wn, input, K, SNRvalues)), 'm'); %Pisarenko
    plot(SNRvalues, log10(errorOnLocations(5, wn, input, K, SNRvalues)), 'color', [1 0.687 0.387]); %Music
    plot(SNRvalues, log10(errorOnLocations(6, wn, input, K, SNRvalues)), 'b'); %Esprit
    plot(SNRvalues, log10(errorOnLocations(7, wn, input, K, SNRvalues)), 'c'); %MatrixPencil
    %%%%%%%% Cramer Rao bound %%%%%%%%%
    
    %Construction of Jacobian Matrix
    F=zeros(N,K);
    for h=1:N
        for j=1:K
            F(h,j)=an(j)*(-1i)*h*exp(-1i*wn(j)*h);
        end
    end
    %Construction of the bound
    inputAmplitude = rms(input);
    ecart_noise = inputAmplitude./10.^(SNRvalues./20);
    cramer_coef = trace(inv(F'*F));
    plot(SNRvalues, log10(cramer_coef.*(ecart_noise.^2)), 'k:');
    
    
    %%%%%%%%% Legends and description %%%%%%%%%%%
    legend('Prony', 'PronyTLS', 'YuleWalker', 'Pisarenko', 'Music', 'Esprit', 'MatrixPencil', 'Cramer Rao');
    descr = {'wn=';wn;'an=';an};
    text(-20,0,descr);
end
