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
%       density: atmospheric density at the given height [N x 1] (g/m^3)

% Variables are persistent so they don't have to be reloaded each time the
% function is called
persistent h  % (m)
persistent rhoMin  % (g/m^3)
persistent rhoMax  % (g/m^3)

% Check if density data has already been allocated
if isempty(h)
    densityData = xlsread("AtmosphericDensity.xlsx");
    h = densityData(:,1);  % (km)
    rhoMin = densityData(:,2);  % (g/km^3)
    rhoMax = densityData(:,3);  % (g/km^3)
    
    % Convert units
    h = h .* 1000;  % (m)
    rhoMin = rhoMin ./ (1000^3);  % (g/m^3)
    rhoMax = rhoMax ./ (1000^3);  % (g/m^3)  
end

% Error checking
if any(r < min(h))
    error("Altitude is below minimum value of 100 km")
end

if any(r > max(h))
   error("Altitude is above maximum value of 1000 km") 
end

% Preallocate
density = zeros(size(r));

% Compute densities
for rIdx = 1:numel(r)
    % 1) Find lower index
    deltas = r(rIdx)-h; % Differences between radius and tabulated heights
    posDeltasIdx = deltas >= 0; % Use > 0 to get lower bound
    [hi, i] = min(deltas(posDeltasIdx)); % Find the smallest difference

    % 2) Calulate scale heights
    Hm = (h(i)-h(i+1)) ./ ( log(rhoMin(i+1)/rhoMin(i)) );
    HM = (h(i)-h(i+1)) ./ ( log(rhoMax(i+1)/rhoMax(i)) );
    
    % 3) Calculate minimum and maximum densities
    densityMin = rhoMin(i)*exp((h(i)-r(rIdx))/Hm);
    densityMax = rhoMax(i)*exp((h(i)-r(rIdx))/HM);
    
    % Compute the density. Ignore diurnal effects.
    density(rIdx) = densityMin;
end