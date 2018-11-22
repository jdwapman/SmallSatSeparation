function [r, w, theta] = trajectory(u)
%trajectory Computes the trajectory of the satellites over time
%   Computes the trajectory of the N satellites over time given the inputs
%   u.
% 
%     Inputs:
%         u: Vector of satellite area commands [N x T] (m^2)
%     Outputs:
%         r: Matrix of satellite altitudes [N x T] (m)
%         w: Matrix of satellite angular velocities [N x T] (rad/sec)
%         theta: Matrix of satellite angular positions [N x T] (rad)
%     globals:
%         N: Number of satellites [1 x 1]
%         T: Number of time steps [1 x 1]
%         r0: Initial radius of the satellites [N x 1] (m)
%         w0: Initial angular velocity of the satellites [N x 1] (rad/sec)
%         theta0: Initial angular position of the satellites [N x 1] (rad)
%         dt: Elapsed time between input commands/calculations [1 x 1] (sec)


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

