%% Setup
clc
clear
load("PhysicalConstants.mat")

%% Test vrel.m

% One dimensional
r = 400E3 + re;  % radius (m)
w = 0.1;  % (rad/sec)
v = vrel(r,w);

handCalc1 = 677100.0;  % (m/sec)
assert(handCalc1 == v)

% Two dimensional
r = [400E3+re; 500E3+re];
w = [0.1; 0.2];
v = vrel(r,w);

handCalc2 = 1374200.0;
handCalc = [handCalc1; handCalc2];
assert(all(v == handCalc))