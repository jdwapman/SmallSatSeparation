function [rNext, wNext, thetaNext] = nextState(r, w, theta, u)
%nextState Calculates the next state of the constellations
%   Calculates the next state of the constellation given the current r, w,
%   and theta vectors as well as the vector of inputs u
%
%   Inputs:
%         r: Current altitude vectors of the satellites [N x 1] (m)
%         w: Current angular velocities of the satellites [N x 1] (rad/sec)
%         theta: Current angular positions of the satellites [N x 1] (rad)
%         u: Satellite area commands [N x 1] (m^2)
%     Outputs:
%         rNext: Next altitude vector of the satellites [N x 1] (m)
%         wNext: Next velocities of the satellites [N x 1] (rad/sec)
%         thetaNext: Next angular positions of the satellites [N x 1] (rad)
%     Globals:
%         dt: Time step between u inputs [1 x 1] (sec)

rNext = r + dt.*Sr(r, w).*u;
wNext = w + dt.*Somega(r,w).*u;
thetaNext = theta + dt.*w + 1/2*dt^2.*Somega(r,w).*u;

end

