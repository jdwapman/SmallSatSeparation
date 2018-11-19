%% Small Satellite Optimization
% This project 1) Replicates the linear-programming-based optimization done
% by Sin et al. and 2) modifies the optimization problem to include nonlinear
% constraints with the goal of evaluating the performance and computational
% complexity of the two optimization methods.
% 
% E. Sin, M. Arcak, and A. Packard, �Small Satellite Constellation Separation
% using Linear Programming based Differential Drag Commands,�
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
N = 105;

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
w0(:) = sqrt(ue./r0.^3);

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
OptimizeLinear()

%% Recreate using nonlinear optimization
OptimizeNonlinear()
