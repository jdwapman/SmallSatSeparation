function [] = OptimizeLinear()
%OptimizeLinear Performs linear optimization of the constellation
%separation problem
%   Performs linear optimization of the constellation separation problem

% Import global variables
global N;  % Number of satellites
global T;  % Number of time steps
global Amin;

% 1) Precompute reference trajectory
u = repmat(Amin, N, T);
[r, w, theta] = trajectory(u);

end

