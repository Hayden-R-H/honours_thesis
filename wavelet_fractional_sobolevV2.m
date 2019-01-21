% @title: Wavelet Threshold for Mix Norms
% @author: Hayden Hohns
% @date: 13/10/16
% @brief: The wavelet transform can be used to 
% detect the number of eigenvalues in the isolated 
% spectrum of an approximation to the transfer and
% composition operators, e.g., taken from Ulam's
% method or similar.

clc;
clear all;

numBins = 1000;
ptsPerBin = 1000;

numev = 30; % for spectral picture
P = ulam(numBins, ptsPerBin);

figure;
spy(P'); % approximate plot of T
xlabel('x_n');
ylabel('Tx_n = x_{n+1}');

[evecR, evalR] = eigs(P, numev, 'LM');
[evecL, evalL] = eigs(transpose(P), numev, 'LM');
[evecR, evalR] = sortem(evecR, evalR);
[evecL, evalL] = sortem(evecR, evalL);

figure;
plot(real(evalR), imag(evalR), 'o', 'MarkerFaceColor', 'b');
xlabel('Re \lambda');
ylabel('Im \lambda');
set(gca, 'fontsize', 18);

% Plot first five eigenvectors
figure;
for i = 1 : 5
    % Plot Koopman (right) eigenvectors
    subplot(2, 5, i);
    plot(real(evecR(:, i)));
    set(gca, 'xtick', []);
    set(gca, 'ytick', []);
    % Plot Frobenius-Perron (left) eigenvectors
    subplot(2, 5, i + 5);
    plot(real(evecL(:, i)));
    set(gca, 'xtick', []);
    set(gca, 'ytick', []);
end

q = 2; % L^q norm
s = 0.4; % smoothing parameter
N = numBins;

maxScale = 600;
scales = 1 : maxScale;

% list of wavelets and markers for final plot
wavelets = {'haar', 'db2', 'db8', 'meyr', ...
    'mexh', 'morl'};
markers = {'*', 'o', 'x', '+', '<', '>'};

sobNormFT = zeros(numev, 1); % Right
sobNormWT = zeros(numev, 1); % Left
figure; 

for k = 1 : length(wavelets)
    
    sobNormFT = zeros(numev, 1); % Left
    waveletName = wavelets{k};
    marker = markers{k};
    
    for i = 1 : numev
        
        l = evecL(:, i); % change for left or right
        %%% Fractional Sobolev norm (Fourier)
        yL = fftshift(fft(l)); % Left
        xi = transpose(- N / 2 : N / 2 - 1);
        zL = yL .* (1 + xi .^ 2) .^ (-s / 2);
        % remove NaN entry
        zL(N / 2 + 1) = yL(N / 2 + 1); 
        sobNormFT(i) = norm(ifft(ifftshift(zL)), q);

        %%% Fractional Sobolev norm via Wavelet Transform
        coeffsL = cwt(l, scales, waveletName); % Left
        energyL = abs(coeffsL) .^ 2;
        temp = sum(energyL, 2);
        
        for j = 1 : length(temp)
            
            sobNormWT(i) = sobNormWT(i) ...
                + temp(j) * (2 ^ (2 * j * s));
            
        end

    end
    
    sobNormFT = sobNormFT / max(sobNormFT);
    sobNormWT = sqrt(sobNormWT);
    sobNormWT = sobNormWT / max(sobNormWT);
    
    %%% Plotting
    hold on;
    plot(sobNormWT, marker, 'MarkerSize', 12);
    set(gca, 'fontsize', 18);
    
    sobNormWT = zeros(numev, 1); % Left
end

%%% Final Plot
hold on;
plot(sobNormFT, '^', 'MarkerSize', 12);
ylabel('Fractional Sobolev Norm');
hh = legend('haar', 'db2', 'db8', 'meyr', 'mexh', ...
    'morl', 'fft');
set(hh, 'Fontsize', 18);