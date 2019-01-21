% @title: 2-D Ulam's method for a variation on the
% double-gyre.
% @author: Hayden Reece Hohns
% @date: 16/11/15
% @brief: An implementation of Ulam's method on a 
% variation of the double-gyre flow using GAIO.

clear all;

% Parameters for the double-gyre flow.
A = 0.25; 
delta = 0.25; 
omega = 2 * pi; 
tStart = 0;

% Double-Gyre flow
v = @(x, T) [-pi * A * x(:, 1) ...
    .* sin(pi * (delta * sin(omega * T) * x(:, 1) .^ 2 ...
    + (1 - 2 * delta * sin(omega * T)) * x(:, 1))) ...
    .* cos(pi * x(:, 2)) ...
    pi * A * cos(pi *(delta * sin(omega * T) ...
    * x(:, 1) .^2 + (1 - 2 * delta * sin(omega * T)) ...
    * x(:, 1))) .* sin(pi * x(:, 2)) .* ...
    (delta * sin(omega * T) * 2 * x(:, 1) ...
    + (1 - 2 * delta * sin(omega * T)))];

n = 100; % number of steps
h = 0.01; % step size
f = @(x) rk4t(v, x, h, n, tStart); % 4th-order Runge-Kutta

nBox = 50; 
x = 2 .* linspace(-1 + 1 / (2 * nBox), ...
    1 - 1 / (2 * nBox), nBox)';
y = (1 / 2) .* x;
% Create grid with 2-cubes for partitions
[XX, YY] = meshgrid(x, y); 
X = [XX(:), YY(:)]; % nBox by nBox set of vectors
c = [1 0.5]; % centre of box
r = [1 0.5]; % radius of box
t = Tree(c, r); % generate empty tree structure
sd = 8; % subdivisions
depth = 15; % Depth of tree.
% Volume of each box is 2^{-depth}.

for i = 1:depth, t.set_flags('all', sd);
    t.subdivide; % make subdivisions in tree structure
end
boxCollection = t.boxes(-1)'; % box collection

% Compute/create transition matrix, P has dimension 2^{13}.
P = tpmatrix(t, f, X, depth, '0')'; % Given in GAIO

% Compute left eigenvectors
[U_l, V_l] = eigs(P', 3);
% Compute (right) eigenvectors
[U_r, V_r] = eigs(P, 3);

figure; % First left eigenvector of U.
show2plus(boxCollection, U_l(:, 1))
title('First Left Eigenvector'); colorbar;

figure; % Second left eigenvector of U.
show2plus(boxCollection, U_l(:, 2))
title('Second Left Eigenvector'); colorbar;