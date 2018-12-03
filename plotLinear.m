function [] = plotLinear(uOptReshape, rMax)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


%% Import Global Variables
global N
global T
global D

%% Plot Angles

% Find actual trajectory and final states
[rAct, wAct, thetaAct] = trajectory(uOptReshape);
t = 1:1:T+1;
rActT = rAct(:,end);  % Extract the final state
wActT = wAct(:,end);
thetaActT = thetaAct(:,end);
figure
polarplot(thetaActT, rActT, '.-')
title("Actual Final State")

%% Plot inputs

figure
t = 1:1:T;
plot(t, uOptReshape);
xlabel("Time (Days)")
ylabel("Area (m^2)")
title("Satellite Area vs Time")

%% Plot radius trajectory

figure
t = 1:1:T+1;
plot(t, rAct/1000)
xlabel("Time (Days)")
ylabel("Altitude (km)")
title("Altitude vs Time")

figure
plot(t, wAct)
xlabel("Time (Days)")
ylabel("Angular Velocity (rad/sec)")
title("Angular Velocity vs Time")

%% Plot angle separation
figure
t = 1:1:T+1;
spacing = rad2deg(D*thetaAct);
plot(t, spacing(1:end-1, :));
xlabel("Time (Days)")
ylabel("Spacing Between Satellites (deg)")
title("Satellite Spacing vs Time")
end

