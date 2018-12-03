function [uOptReshape, rMax] = OptimizeLinear(T, rStart, wStart, thetaStart)
%OptimizeLinear Performs linear optimization of the constellation
%separation problem
%   Performs linear optimization of the constellation separation problem

%% Import global variables

% IMPORTANT NOTE: DO NOT IMPORT THE GLOBAL T, r0, w0, or theta0 variables.
% Since this optimization function is used for the closed-loop
% model-predictive control version, the starting variables may change from
% iteration to iteration. Instead, initial conditions should be passed to
% the function as rStart, wStart, and thetaStart. Additionally, T (the time
% horizon) will change from iteration to iteration when performing model
% predictive control

global Amin;
global Amax;
global N;
global dt;
global D;
global delta_des;
global epsTheta;
global epsOmega;
global re;

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

b = [rStart; ...
    epsTheta * ones(N,1) - D*(thetaStart + dt*T*wStart) + delta_des; ...
    epsTheta * ones(N,1) + D*(thetaStart + dt*T*wStart) - delta_des; ...  
    epsOmega * ones(N,1) - D*wStart; ...
    epsOmega * ones(N,1) + D*wStart; ...
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
options = optimoptions('fmincon', 'Display', 'iter', 'MaxFunctionEvaluations', 100*N*T, 'UseParallel', true);
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
assert(all(rMax < rStart))  % Check that the radius is less than starting

end

