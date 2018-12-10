function [uOptReshape, rMax] = OptimizeLinear()
%OptimizeLinear Performs linear optimization of the constellation
%separation problem
%   Performs linear optimization of the constellation separation problem
%
%   Inputs:
%       None
%   Outputs:
%       uOptReshape: Satellite area commands for N satellites across T time
%           steps [N x T] (m^2)
%       rMax: The maximized minimum satellite radius (m)
%   Globals:
%       Almost all globals are imported from the SmallSatOptimization
%       function. For definitions, refer to that file. In general it is not
%       ideal to use global functions, but due to the nature of Matlab's
%       nonlinear constraints functions in fmincon (the function must be
%       called many times and arguments cannot be passed to it), all
%       physical and simulation parameters must be in the global
%       workspace. Writing to and reading from a file is too slow.

%% Import global variables

global Amin;
global Amax;
global N;
global dt;
global D;
global delta_des;
global epsTheta;
global epsOmega;
global re;
global T;
global r0;
global w0;
global theta0;

%% Precompute reference trajectory
uMin = repmat(Amin, N, T);
[rRef, wRef, thetaRef] = trajectory(uMin);

% Precompute Sr and Somega matrices
Sr_ref = [];
Somega_ref = [];
Salpha_ref = [];
alpha = ((T-1/2)*ones(1,T) - (0:1:T-1));
for n = 1:1:N
    Sr_i = Sr(rRef(n,1:T), wRef(n,1:T))';  % Must convert to row vector
    Sr_ref = blkdiag(Sr_ref, Sr_i);  % Append matrix to diagonal
    
    Somega_i = Somega(rRef(n,1:T), wRef(n,1:T))';  % Must convert to row vector
    Somega_ref = blkdiag(Somega_ref, Somega_i);
    
    Salpha_i = alpha .* Somega_i;
    Salpha_ref = blkdiag(Salpha_ref, Salpha_i);
end


%% Create  matrices for linear programming
A = [-dt * Sr_ref, -ones(N,1); ...
    dt^2 * D * Salpha_ref, zeros(N,1); ...
    -dt^2 * D * Salpha_ref, zeros(N,1); ...
    dt * D * Somega_ref, zeros(N,1); ...
    -dt * D * Somega_ref, zeros(N,1); ...
    eye(N*T, N*T), zeros(N*T, 1); ...
    -eye(N*T, N*T), zeros(N*T, 1)];

b = [r0; ...
    epsTheta * ones(N,1) - D*(theta0 + dt*T*w0) + delta_des; ...
    epsTheta * ones(N,1) + D*(theta0 + dt*T*w0) - delta_des; ...  
    epsOmega * ones(N,1) - D*w0; ...
    epsOmega * ones(N,1) + D*w0; ...
    Amax * ones(N*T,1); ...
    -Amin * ones(N*T,1)];

%% Create cost function
f = [zeros(1,N*T), 1];
x0 = [repmat(Amin, N*T, 1); -(100e3+re)];  % Initial guess

costfun = @(x)f*x;

Aeq = [];
beq = [];
lb = [];
ub = [];
nonlcon = [];

%% Optimize

% Perform optimization. NOTE: for final results, disable parallelization
options = optimoptions('fmincon', 'Display', 'iter', 'MaxFunctionEvaluations', 1000*N*T, 'UseParallel', false);
result = fmincon(costfun, x0, A, b, Aeq, beq, lb, ub, nonlcon, options);

uOpt = result(1:end-1);  % Extract area commands
thresh = result(end);  % Extract radius

% IMPORTANT NOTE: We originally reshaped to N*T and neglected to perform
% the transpose operation, which led to incorrect results. The "results"
% vector from the optimization function is [N*T, 1] where the form is in
% the following format in a vector
% N = 1 [Time commands T x 1]
% N = 2 [Time commands T x 1]
% etc...
% Reshaping must be first be done to T x N, then the transpose must be
% taken so that the uOptReshape matrix is an [N x T] matrix where rows are
% satellites and columns are commands at a given time step

uOptReshape = reshape(uOpt, T, N);
uOptReshape = uOptReshape';

rMax = -thresh;  % Determines maximized radius
assert(all(rMax < r0))  % Check that the radius is less than starting

end

