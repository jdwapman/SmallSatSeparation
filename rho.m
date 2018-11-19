function [density] = rho(r)
%rho Returns the atmospheric density at the altitude r
%   Uses a simplified version of the Harris-Priester model, as described by
%   Montenbruck, that does not include diurnal effects.
%
%   O. Montenbruck and E. Gill, Satellite Orbits. New York: Springer,
%   2005, ch. 3.
%
%   Currently, we do not have access to this book so this function is
%   a placeholder.
%
%   Inputs:
%       r: Satellite altitude [N x 1] (m)
%   Outputs:
%       density: atmospheric density at the given height [N x 1] (unknown)


% TODO: Get actual equation
density = ones(size(r));

end

