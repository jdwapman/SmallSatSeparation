function [r, w, theta] = trajectoryLinear(u, Sr_ref, Somega_ref, Salpha_ref)
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
global N;
global T;
global r0;
global w0;
global theta0;
global dt;

% Preallocate trajectory outputs
r = zeros(N, T);
w = zeros(N, T);
theta = zeros(N, T);

uflat = reshape(u, N*T, 1);

% Iterate through all other time steps. Functions from Sin et al. paper
for t = 1:1:T
    r(:,t) = r0 + dt*Sr_ref(:,1:N*t)*uflat(1:N*t);
    w(:,t) = w0 + dt*Somega_ref(:,1:N*t)*uflat(1:N*t);
    theta(:,t) = theta0+dt*t*w0+dt^2*Salpha_ref(:,1:N*t)*uflat(1:N*t);
end

r = [r0 r];
w = [w0 w];
theta = [theta0 theta];

end

