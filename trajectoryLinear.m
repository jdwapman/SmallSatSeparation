function [r, w, theta] = trajectoryLinear(u, rRef, wRef, thetaRef)
%trajectory Computes the linearized trajectory of the satellites over time
%   Computes the linearized trajectory of the N satellites over time
%   given the inputs u and the reference trajectory.
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
% Note: these are [N x T+1] to include initial state
uRef = repmat(Amin, N, T);
[rRef, wRef, thetaRef] = trajectory(uRef);

% Preallocate linear trajectory outputs
r = zeros(N, T);
w = zeros(N, T);
theta = zeros(N, T);


% Iterate through all other time steps. Functions from Sin et al. paper
for t = 1:1:T

    % Compute for r
    rSum = 0;
    for k = 0:1:t-1
       rSum = rSum + Sr(rRef(:,k+1), wRef(:,k+1)).* u(:,k+1);
    end

    r(:,t) = r0 + dt*rSum;
    
    % Compute for w
    wSum = 0;
    for k = 0:1:t-1
       wSum = wSum + Somega(rRef(:,k+1), wRef(:,k+1)).* u(:,k+1);
    end

    w(:,t) = w0 + dt*wSum;
   
    
end

% NOTE: NOT ACCURATE
for t = 1:1:T
    % Compute for theta
    wThetaSum = 0;
    thetaSum = 0;
    for k = 0:1:t-1
       wThetaSum = wThetaSum + w(:,k+1);
       thetaSum = thetaSum + Somega(rRef(:,k+1), wRef(:,k+1)).* u(:,k+1);
    end

    theta(:,t) = theta0 + dt*wThetaSum + 0.5*(dt^2)*thetaSum;
end

% Append initial position
r = [r0 r];
w = [w0 w];
theta = [theta0 theta];

end

