clc;
clear all;

numBins = 1000;
ptsPerBin = 1000;

numev = 100; % for spectral picture
P = ulam(numBins, ptsPerBin);

figure;
spy(P'); % approximate plot of T for illustration purposes
title('Map');
xlabel('x_n');
ylabel('Tx_n = x_{n+1}');

[evecR, evalR] = eigs(P, numev, 'LM');
[evecL, evalL] = eigs(transpose(P), numev, 'LM');
[evecR, evalR] = sortem(evecR, evalR);
[evecL, evalL] = sortem(evecR, evalL);
figure;
plot(real(evalR), imag(evalR), 'o');
title('Spectrum of the Koopman Operator');
xlabel('Re \lambda');
ylabel('Im \lambda');

% Plot first five eigenvectors
figure;
for i = 1 : 5
    % Plot Koopman (right) eigenvectors
    subplot(2, 5, i);
    plot(real(evecR(:, i)));
    % Plot Frobenius-Perron (left) eigenvectors
    subplot(2, 5, i + 5);
    plot(real(evecL(:, i)));
end

q = 2; % L^q norm
s = 0.4; % smoothing parameter
N = numBins;

maxScale = 700;
scales = 1 : maxScale;

% change for different wavelet choice below
waveletName = 'haar';

sobNormFT = zeros(numev, 1); % Right
sobNormWT = zeros(numev, 1); % Right
sobNormFTL = sobNormFT; % Left
sobNormWTL = sobNormWT; % Left

for i = 1 : numev
    w = evecR(:, i); % i-th right eigenvector
    l = evecL(:, i); % i-th left eigenvector
    
    %%% Fractional Sobolev norm via Fourier Transform (RIGHT)
    y = fftshift(fft(w)); % Right
    xi = transpose(- N / 2 : N / 2 - 1);
    z = y .* (1 + xi .^ 2) .^ (-s / 2);
    z(N / 2 + 1) = y(N / 2 + 1); % remove NaN entry in z
    sobNormFT(i) = norm(ifft(ifftshift(z)), q);
    
    %%% Fractional Sobolev norm via Wavelet Transform (RIGHT)
    coeffsR = cwt(w, scales, waveletName); % Right
    energyR = abs(coeffsR) .^ 2;
    temp = sum(energyR, 2);
    for j = 1 : length(temp)
        sobNormWT(i) = sobNormWT(i) ...
            + temp(j) * (2 ^ (2 * j * s));
    end
    
    %%% Fractional Sobolev norm via Fourier Transform (LEFT)
    yL = fftshift(fft(l)); % Left
    zL = yL .* (1 + xi .^ 2) .^ (-s / 2);
    zL(N / 2 + 1) = yL(N / 2 + 1); % remove NaN entry in z
    sobNormFTL(i) = norm(ifft(ifftshift(zL)), q);
    
    %%% Fractional Sobolev norm via Wavelet Transform (LEFT)
    coeffsL = cwt(l, scales, waveletName); % Left
    energyL = abs(coeffsL) .^ 2;
    temp = sum(energyL, 2);
    for j = 1 : length(temp)
        sobNormWTL(i) = sobNormWTL(i) ...
            + temp(j) * (2 ^ (2 * j * s));
    end
    
end

sobNormFT = sobNormFT / max(sobNormFT);
sobNormFTL = sobNormFTL / max(sobNormFTL);


sobNormWT = sqrt(sobNormWT);
sobNormWTL = sqrt(sobNormWTL);

sobNormWT = sobNormWT / max(sobNormWT);
sobNormWTL = sobNormWTL / max(sobNormWTL);

%%% Plotting
figure; 
plot(sobNormWT, '^', 'MarkerSize', 12);
xlabel('Right and Left Eigenvectors');
ylabel('Fractional Sobolev Norm');
set(gca, 'fontsize', 18);
hold on;
plot(sobNormWTL, 'v', 'MarkerSize', 12);
hold on;
plot(sobNormFT, '*', 'MarkerSize', 12);
hold on;
plot(sobNormFTL, 'o', 'MarkerSize', 12);
hh = legend('Right Wavelet', 'Left Wavelet', ...
    'Right Fourier', 'Left Fourier');
set(hh, 'Fontsize', 12);