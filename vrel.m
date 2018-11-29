function [result] = vrel(r,w)
%vrel Returns the angular velocity of a satellite relative to the earth
%   Computes the angular velocity magnitude of a satellite relative
%   to the earth's angular rotation
%
%   Inputs:
%       r: Satellite orbital radius [N x 1] (m)
%       w: Satellite angular velocities [N x 1] (rad/sec)
%   Outputs:
%       result: satellite velocity relative to the earth's angular rotation
%   Globals:
%       we: Angular velocity of the earth (rad/sec)
%       phi: inclination of the satellite's orbit.
%               90 deg: near-polar
%               0 deg: equitorial
%               180 deg: Retrograde equitorial

% Need to declare needed globals within the function
global we
global phi

result = r.*(w-we*cos(phi));

end

