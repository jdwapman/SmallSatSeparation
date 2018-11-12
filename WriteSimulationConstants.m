function [] = WriteSimulationConstants(filename)
%WriteSimulationConstants Writes simulation constants to constants.mat
%   The nonlinear optimization functions all need access to the various
%   simulation parameters. This function writes these parameters to
%   constants.mat so they can be easily read from within any function
%   simply by loading the mat file into the workspace

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
    
    save(filename)
end

