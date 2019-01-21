% @title: Wavelet threshold for the quadruple-gyre
% @author: Hayden Hohns
% @date: 31/10/16
% @brief: An application of the wavelet fractional
% Sobolev norm algorithm for detecting the number 
% of elements in  the isolated spectrum. This 
% algorithm requires the left and right eigenvectors
% of an Ulam or DMD matrix.

clc;
clear all;

tStart = 0;

% vector field parameters
delta = 0.1;
omega = pi / 5;
A = 0.1;
B = 1;

% vector field
ff = @(t, x) delta * sin(omega * t) .* x .^ 2 ...
    + (1 - delta * sin(omega * t)) .* x;
dff = @(t, x) delta * sin(omega * t) .* x + 1 ...
    - B * delta * sin(omega * t);
g = @(t, x, y) sin(A * pi * ff(t, x)) ...
    .* cos(B * pi * ff(t, y)) .* dff(t, y);
v = @(x, T) [-g(T, x(:, 1), x(:, 2)) ...
    g(T, x(:, 2), x(:, 1))];

% n is the number of time steps, h is the step size.
n = 100;
h = 0.01;
% 4th-order Runge-Kutta.
f = @(x) rk4t(v, x, h, n, tStart); 

nBox = 50; 
x = 2 .* linspace(-1 + 1 / (2 * nBox), ...
    1 - 1 / (2 * nBox), nBox)';
y = (1 / 2) .* x;
% Create grid with 2-cubes for partitions.
[XX, YY] = meshgrid(x, y); 
X = [XX(:), YY(:)]; % nBox by nBox set of vectors.
c = [1 0.5]; % Centre of box.
r = [1 0.5]; % Radius of box.
t = Tree(c, r); % Generate empty tree structure.
sd = 8; % Subdivisions.
depth = 12; % Depth of tree.
% Volume of each box is 2^{-depth}.

for i = 1:depth, t.set_flags('all', sd);
    % Make subdivisions in tree structure.
    t.subdivide; 
end
boxCollection = t.boxes(-1)';

% Compute transition matrix, P has dimension 2^{12}.
P = tpmatrix(t, f, X, depth, '0')'; % Given in GAIO

K = 20;
% Compute left eigenvectors
[U_l, V_l] = eigs(P', K);
% Compute (right) eigenvectors
[U_r, V_r] = eigs(P, K);

% sort eigenvalues and eigenvectors
[U_l, V_l] = sortem(U_l, V_l);
[U_r, V_r] = sortem(U_r, V_r);

figure; % First left eigenvector of U.
show2plus(boxCollection, U_l(:, 1))
title('First Left Eigenvector'); 
colorbar;

figure; % Second left eigenvector of U.
show2plus(boxCollection, U_l(:, 2))
title('Second Left Eigenvector'); 
colorbar;

figure; % Third left eigenvector of U.
show2plus(boxCollection, U_l(:, 2))
title('Third Left Eigenvector'); 
colorbar;

figure; % Fourth left eigenvector of U.
show2plus(boxCollection, U_l(:, 2))
title('Fourth Left Eigenvector'); 
colorbar;

q = 2; % L^q norm
s = 1; % smoothing parameter

maxScale = 120;
scales = 1 : maxScale;
waveletName = 'sinc';
sobNormWT = zeros(K, 1); % Right
sobNormWTL = sobNormWT; % Left

for i = 1 : K
    
    % i-th right eigenvector
    w = reshape(U_r(:, i), [128, 32]); 
    % i-th left eigenvector
    l = reshape(U_l(:, i), [128, 32]); 
    [M, N] = size(w);
    
    %%% Wavelet Fractional Sobolev Norm
    % Right
    coeffsR = cwtft2(w, 'wavelet', waveletName, ...
        'scales', scales, 'angles', 0); 
    energyR = abs(coeffsR.cfs) .^ 2; % from struct
    temp1 = sum(energyR, 2);
    temp2 = sum(temp1, 3);
    [cols, rows] = size(temp2);
    
    for j = 1 : cols
        for k = 1 : rows 
            sobNormWT(i) = sobNormWT(i) ...
                + temp2(j, k) * 2 ^ (-j * k * s);
        end
    end
    
    % Left
    coeffsL = cwtft2(l, 'wavelet', waveletName, ...
        'scales', scales, 'angles', 0); 
    energyL = abs(coeffsL.cfs) .^ 2; % from struct
    temp1 = sum(energyL, 2);
    temp2 = sum(temp1, 3);
    [cols, rows] = size(temp2);
    
    for j = 1 : cols
        for k = 1 : rows
            sobNormWTL(i) = sobNormWTL(i) ...
                + temp2(j, k) * 2 ^ (-j * k * s);
        end
    end
    
end

sobNormWT = sobNormWT / max(sobNormWT);
sobNormWTL = sobNormWTL / max(sobNormWTL);

%%% Plotting
figure; 
plot(sobNormWT, '*', 'MarkerSize', 12);
xlabel('Right eigenvectors');
ylabel('Fractional Sobolev Norm');

figure; 
plot(sobNormWTL, '*', 'MarkerSize', 12);
xlabel('Left eigenvectors');
ylabel('Fractional Sobolev Norm');