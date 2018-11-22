function [r, w, theta] = trajectory(u)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Import global variables
global N; % Number of satellites
global T; % Number of time steps
global r0;
global w0;
global theta0;
global dt;

% Preallocate trajectory outputs
r = zeros(N, T);
w = zeros(N, T);
theta = zeros(N, T);

% Set column 1 to initial conditions
r(:,1) = r0;
w(:,1) = w0;
theta(:,1) = theta0;

for t = 2:1:T
    r(:,t) = r(:,t-1) + (dt.*Sr(r(:,t-1), w(:,t-1)) .* u(:,t-1));
    w(:,t) = w(:,t-1) + (dt.*Somega(r(:,t-1), w(:,t-1)) .* u(:,t-1));
    theta(:,t) = theta(:,t-1) + dt.*w(:,t-1) + (0.5 * dt^2 .*Somega(r(:,t-1), w(:,t-1)) .* u(:,t-1));
end

end

