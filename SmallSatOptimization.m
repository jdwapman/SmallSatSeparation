%% Small Satellite Optimization
% DESCRIPTIVE TEXT

%% Initialization Steps
clc
clear

%% Set Up Constants

% Earth's gravitational constant. (Wikipedia)
ue = 3.986004418e14;  % m^3s^-2
GM = ue;  % Alternative representation

% Satellite drag coefficient. (Li & Mason)
Cd = 2.2;  

% Surface area exposed to the incident stream. (Li & Mason)
% Amin < A < Amax
Amax = 0.225;  % (m^2)
Amin = 0.0371; % (m^2)

% Satellite mass. (Li & Mason)
m = 4.9;  % (kg)

% Velocity of satellite relative to the atmosphere
vrel = 0;  % (m/sec)

% Atmospheric density at the satellite position (kg/m^3)
rho = 2.5e-12;  

%% Simulation Parameters

% Number of satellites
N = 105; 

% Initial theta
theta_0 = zeros(N,1);

% r_0
% TODO

% omega_0
% TODO

% Circular orbit altitude
h = 475 * 10^3;  % (m)

% Theta tolerance
eps_theta = 0.1;  % Degrees
eps_omega = 1e-18;  % Rad/sec

 