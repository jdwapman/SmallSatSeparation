function [result] = Somega(r,w)
%Somega Computes the Sr term needed for the satellite orbital dynamics
%     Computes the Sr term needed for the satellite orbital dynamics.
%
%     "The values for SR() and SW() describe how the impact of
%     the input changes depending on the current state of the satellite.
%     For example, for any given cross-sectional surface area,
%     a satellite will experience greater atmospheric drag force at
%     lower altitudes than at higher altitudes. This relationship is
%     captured using the Gaussian variation of parameters (VOP)
%     form of the equations of motion. These equations are used to
%     approximate the rates of change of the time-varying elements
%     in the solution for the unperturbed, two-body system due to
%     small perturbing forces. Vallado [15] shows that the average
%     rate of change in the semi-major axis of an orbit and the
%     angular speed of the satellite can be expressed in terms of the
%     atmospheric drag perturbation." from
% 
%     E. Sin, M. Arcak, and A. Packard, “Small Satellite Constellation
%     Separation using Linear Programming based Differential Drag
%     Commands,” arXiv:1710.00104 [cs], Sep. 2017.
% 
%     Inputs:
%         r: Satellite altitude [N x 1] (m)
%         w: Satellite angular velocity [N x 1] (rad/sec)
%     Outputs:
%         result: Somega term
%     Globals:
%         Cd: Satellite drag coefficient [1 x 1]
%         m: Satellite mass [N x 1] or [1 x 1] (kg)
%         ue: Earth's gravitational constant (m^3s^-2)

global Cd
global m
global ue

result = 3/2.*Cd./m.*rho(r).*abs(vrel(r,w)).^2./r;

end

