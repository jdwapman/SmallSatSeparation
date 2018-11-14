%% Small Satellite Optimization
% Project Description

%% Initialization Steps
clc
clear

%% Set Up Constants

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

% Angular velocity of the earth TODO: Get value
global we;
we = 1;  % (rad/sec)

% Atmospheric density at the satellite position. TODO: Get Montenbruck book
global rho;

%% Simulation Parameters

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

%% Recreate using linear programming

% Creates D matrix used for satellite spacing
global D;
D = eye(N,N);
D = D - diag(ones(N-1,1),1);
D(N,1) = -1;

% Creates delta_des
global delta_des;
delta_des = repmat(2*pi/N, N, 1);
delta_des(end) = -2*pi/N*(N-1);  % Replace last value

%% Recreate using nonlinear optimization

