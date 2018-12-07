function [density] = rho(r)
%rho Returns the atmospheric density at the altitude r
%   Uses a simplified version of the Harris-Priester model, as described by
%   Montenbruck, that does not include diurnal effects.
%
%   O. Montenbruck and E. Gill, Satellite Orbits. New York: Springer,
%   2005, ch. 3.
%
%   Inputs:
%       r: Satellite orbital radius [N x 1] (m)
%   Outputs:
%       density: atmospheric density at the given height [N x 1] (kg/m^3)

% Variables are persistent so they don't have to be reloaded each time the
% function is called
persistent h  % (m)
persistent rhoMin  % (kg/m^3)
persistent rhoMax  % (kg/m^3)

global re;  % Radius of the earth (m)

% Check if density data has already been allocated
if isempty(h)
    densityData = xlsread("AtmosphericDensity.xlsx");
    hExcel = densityData(:,1);  % (km)
    rhoMinExcel = densityData(:,2);  % (g/km^3)
    rhoMaxExcel = densityData(:,3);  % (g/km^3)
    
    % Convert units
    h = hExcel .* 1000;  % (m)
    rhoMin = rhoMinExcel ./ (1000^3) ./ 1000;  % (kg/m^3)
    rhoMax = rhoMaxExcel ./ (1000^3) ./ 1000;  % (kg/m^3)  
end

altitude = r - re;

if any(altitude < 0)
   error("Altitude < 0") 
end

% Error checking
if any(altitude < 100e3) % 100 km
    warning("Altitude is below minimum value of 100 km")
end

if any(altitude > 1000e3) % 1000 km
   warning("Altitude is above maximum value of 1000 km") 
end


% Preallocate
density = zeros(size(altitude));

% Compute densities
for rIdx = 1:numel(altitude)
    % 1) Find lower index
    deltas = altitude(rIdx)-h; % Differences between radius and tabulated heights
    posDeltasIdx = deltas >= 0; % Use > 0 to get lower bound
    [hi, i] = min(deltas(posDeltasIdx)); % Find the smallest difference

    % 2) Calulate scale heights
    Hm = (h(i)-h(i+1)) ./ ( log(rhoMin(i+1)/rhoMin(i)) );
    HM = (h(i)-h(i+1)) ./ ( log(rhoMax(i+1)/rhoMax(i)) );
    
    % 3) Calculate minimum and maximum densities
    densityMin = rhoMin(i)*exp((h(i)-altitude(rIdx))/Hm);
    densityMax = rhoMax(i)*exp((h(i)-altitude(rIdx))/HM);
    
    % 4) Return the density. Ignore diurnal effects.
    density(rIdx) = densityMin;
end

% density(:) = 2.5E-12; % Constant density value from Li and Mason. Useful for debugging