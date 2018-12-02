function [r, w, theta] = trajectoryLinear(u)
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
global Amin;

% 1) Precompute reference trajectory
uRef = repmat(Amin, N, T);
[rRef, wRef, thetaRef] = trajectory(uRef);  % Note: these are [N x T+1] to include initial state

% Preallocate linear trajectory outputs
r = zeros(N, T);
w = zeros(N, T);
theta = zeros(N, T);

% Iterate through all other time steps. Functions from Sin et al. paper
for t = 1:1:T

    rSum = 0;
    for sumT = 0:1:t-1
       rSum = rSum + Sr(rRef(:,sumT+1), wRef(:,sumT+1)).*u(:,sumT+1);
    end
    
    r(:,t) = r0 + dt*rSum;
    
end

r = [r0 r];

w = 0;
theta = 0;

end

