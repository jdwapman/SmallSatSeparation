function [u,rMax] = OptimizeMPC(method)
%OptimizeLinearMPC Summary of this function goes here
%   Detailed explanation goes here


%% Import global variables

global N;
global T;
global r0;
global w0;
global theta0;

%% Preallocate matrices
r = zeros(N, T+1);
w = zeros(N, T+1);
theta = zeros(N, T+1);
u = zeros(N, T);

%% Set first time step to initial conditions
r(:,1) = r0;
w(:,1) = w0;
theta(:,1) = theta0;

% Perform MPC optimization

fullHorizon = T+1;
for t = 1:1:fullHorizon-1
    
    % Execute optimization problem using the last state as the initial
    % conditions. Also, shrink the T horizon
    T = fullHorizon - t;
    
    disp(T)
    
    r0 = r(:,t);
    w0 = w(:,t);
    theta0 = theta(:,t);
    
    if strcmp(method, 'linear')
        [commands, rMax] = OptimizeLinear();
    elseif strcmp(method, 'nonlinear')
        [commands, rMax] = OptimizeNonlinear();
    else
        error("Invalid Optimization Method")
    end
    
    [rCalc, wCalc, thetaCalc] = trajectory(commands);
    
    % Position 1 stores initial conditions. Take position 2 for the next
    % step
    r(:,t+1) = rCalc(:,2); 
    w(:,t+1) = wCalc(:,2);
    theta(:,t+1) = thetaCalc(:,2); 
    
    u(:,t) = commands(:,1);
end

% Restore original values
T = fullHorizon - 1;
r0 = r(:,1);
w0 = w(:,1);
theta0 = theta(:,1);

rMax = min(r(:,end));

end

