%% Small Satellite Optimization
% This project 1) Replicates the linear-programming-based optimization done
% by Sin et al. and 2) modifies the optimization problem to include nonlinear
% constraints with the goal of evaluating the performance and computational
% complexity of the two optimization methods.
% 
% E. Sin, M. Arcak, and A. Packard, “Small Satellite Constellation Separation
% using Linear Programming based Differential Drag Commands,”
% arXiv:1710.00104 [cs], Sep. 2017.

%% Initialization Steps
clc
clear

%% Set Up Physical Constants

% Note, must use global variables so they are available to the functions in
% the nonlinear constraints. Keeping these variables in memory is faster
% than saving to a file

% Earth's gravitational constant. (Wikipedia)
global ue;
ue = 3.986004418e14;  % m^3s^-2

% Satellite drag coefficient. (Li & Mason)
global Cd;
Cd = 2.2;  

% Surface area exposed to the incident stream. (Li & Mason)
% Amin < A < Amax
global Amin Amax
Amin = 0.0371; % (m^2)
Amax = 0.225;  % (m^2)

% Satellite mass. (Li & Mason)
global m;
m = 4.9;  % (kg)

% Angular velocity of the earth. (Wikipedia)
global we;
we = 7.2921150e-5;  % (rad/sec)

% Inclination of orbit
global phi;
phi = deg2rad(90);  % 90 degrees for near-polar orbit


%% Simulation Parameters
% Specified in Sin et al. paper

% Number of satellites in the simulation
global N;
N = 105; % 105;

% Number of time steps/input commands
global T;
T = 71;  % (Days)

% Initial position
global theta0;
theta0 = zeros(N,1);  % Preallocate for speed
theta0(:) = 0;      % Set to initial conditions

% Initial altitude
global r0;
r0 = zeros(N,1);  % Preallocate for speed
r0(:) = 475000;   % (m)

% Initial velocity
global w0;
w0 = zeros(N,1);
w0(:) = sqrt(ue./(r0.^3));

% Update time
global dt;
dtDays = 1;
dt = dtDays * 24 * 60 * 60;

% Position tolerance
global epsTheta;
epsTheta = 0.1;  % (rad)

% Velocity tolerance
global epsOmega;
epsOmega = 1e-18;  % (rad/sec)

%% Variables for both programming methods

% Creates D matrix used for satellite spacing
global D;
D = eye(N,N);
D = D - diag(ones(N-1,1),1);
D(N,1) = -1;

% Creates delta_des
global delta_des;
delta_des = repmat(2*pi/N, N, 1);
delta_des(end) = -2*pi/N*(N-1);  % Replace last value

%% Recreate using linear programming

% OptimizeLinear()

% 1) Precompute reference trajectory
u = repmat(Amin, N, T);
[rRef, wRef, thetaRef] = trajectory(u);

% 2) Precompute Sr and Somega matrices
Sr_ref = [];
Somega_ref = [];
Salpha_ref = [];
alpha = ((T-1/2)*ones(1,T) - (0:1:T-1));
for n = 1:1:N
    Sr_i = Sr(rRef(n,:), wRef(n,:))';  % Must convert to row vector
    Sr_ref = blkdiag(Sr_ref, Sr_i);
    
    Somega_i = Somega(rRef(n,:), wRef(n,:))';  % Must convert to row vector
    Somega_ref = blkdiag(Somega_ref, Somega_i);
    
    Salpha_i = alpha .* Somega_i;
    Salpha_ref = blkdiag(Salpha_ref, Salpha_i);
end


% 4) Create A matrix for linear programming
A = [-dt*Sr_ref, -ones(N,1); ...
    dt^2*D*Salpha_ref, zeros(N,1); ...
    -dt^2*D*Salpha_ref, zeros(N,1); ...
    dt*D*Somega_ref, zeros(N,1); ...
    -dt*D*Somega_ref, zeros(N,1); ...
    eye(N*T, N*T), zeros(N*T, 1); ...
    -eye(N*T, N*T), zeros(N*T, 1)];

b = [r0; ...
    epsTheta * ones(N,1) - D*(theta0 + dt*T*w0)+delta_des; ...
    epsTheta * ones(N,1) + D*(theta0 + dt*T*w0)-delta_des; ...  
    epsOmega * ones(N,1) - D*w0; ...
    epsOmega * ones(N,1) + D*w0; ...
    Amax * ones(N*T, 1); ...
    -Amin * ones(N*T, 1)];

% Create cost function
f = [zeros(1,N*T), 1];
x0 = [repmat(Amin, N*T, 1); 100e3];

costfun = @(x)f*x;

% Optimize
result = fmincon(costfun, x0, A, b);

% Plot reference trajectories
% t = 1:1:T;
% figure
% plot(t,rRef/1000)
% xlabel("Time (Days)")
% ylabel("Altitude (km)")
% title("Altitude vs Time")
% 
% figure
% plot(t, rad2deg(wRef))
% xlabel("Time (Days)")
% ylabel("Angular Velocity (deg)")
% title(" Angular Velocity vs Time")


%% Recreate using nonlinear optimization
OptimizeNonlinear()
