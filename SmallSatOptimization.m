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
clc  % Clear the command line
clear  % Remove all variables from the workspace
close all  % Close all figures

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

% Inclination of orbit.  (Sin)
global phi;
phi = deg2rad(90);  % 90 degrees for near-polar orbit

% Radius of the earth
global re;
re = 6371E3;  % (m)

save("Results/PhysicalConstants.mat")

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
r0(:) = 475000 + re;   % Convert from altitude to radius(m)

% Initial velocity
global w0;
w0 = zeros(N,1);
w0(:) = sqrt(ue./(r0.^3));  % (rad/sec) TODO: Needs to be multiplied by 2pi? Current units are 1/sec

% Update time
global dt;
dtDays = 1;
dt = dtDays * 24 * 60 * 60;  % (sec)

% Position tolerance
global epsTheta;
epsTheta = deg2rad(0.1);  % (rad)

% Velocity tolerance
global epsOmega;
epsOmega = 1e-18;  % (rad/sec)

% Atmospheric density settings
global densitySetting;
densitySetting = "low";

save("Results/SimulationParameters.mat")

%% Variables for all programming methods

% Creates D matrix used for satellite spacing
global D;
D = eye(N,N);
D = D - diag(ones(N-1,1),1);
D(N,1) = -1;  % Bottom-left corner

% Creates delta_des
global delta_des;
delta_des = repmat(2*pi/N, N, 1);
delta_des(end) = -2*pi/N*(N-1);  % Replace last value

%% Choose modes to calculate
linear = 1;
linearMPC = 0;
nonlinear = 0;
nonlinearMPC = 0;

%% Recreate using open-loop linear programming

if linear
    % Start timing
    tic

    [commandsLinearOpenLoop, rLinearOpenLoop] = OptimizeLinear();

    % Stop timing
    tLinearOpenLoop = toc;
    
    
    filename = ['Results/Linear', 'N', string(N), 'T', string(T), '.mat'];
    filename = strjoin(filename, '');
    save(filename);
end

%% Plot linear optimization results

if linear
    plotSatellites(commandsLinearOpenLoop, rLinearOpenLoop);
end

%% Recreate using closed-loop linear programming (model-predictive control)

if linearMPC
    % Start timing
    tic

    [commandsLinearClosedLoop, rLinearClosedLoop] = OptimizeMPC('linear');

    % Stop timing
    tLinearClosedLoop = toc;
    
    filename = ['Results/LinearMPC', 'N', string(N), 'T', string(T), '.mat'];
    filename = strjoin(filename, '');
    save(filename);
end

%% Plot

if linearMPC
    plotSatellites(commandsLinearClosedLoop, rLinearClosedLoop);
end
%% Nonlinear open-loop optimization

if nonlinear
    % Start timing
    tic

    [commandsNonlinearOpenLoop, rNonlinearOpenLoop] = OptimizeNonlinear();

    tNonlinearOpenLoop = toc;
    
    filename = ['Results/NonlinearOpenLoop', 'N', string(N), 'T', string(T), '.mat'];
    filename = strjoin(filename, '');
    save(filename);
end

%% Plot Nonlinear open-loop optimization

if nonlinear
    plotSatellites(commandsNonlinearOpenLoop, rNonlinearOpenLoop)
end

%% Nonlinear closed-loop optimization

if nonlinearMPC
    tic

    [commandsNonlinearClosedLoop, rNonlinearClosedLoop] = OptimizeMPC('nonlinear');

    tNonlinearClosedLoop = toc;
    
    filename = ['Results/NonlinearMPC', 'N', string(N), 'T', string(T), '.mat'];
    filename = strjoin(filename, '');
    save(filename);
end