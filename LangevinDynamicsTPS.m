function [outPath, outZ] = LangevinDynamicsTPS(ZIn, idx, dt, tf, eta, Gamma, kT, h4, h2)
    % This function simulates the Langevin dynamics for a system using two-point sampling (TPS).
    % Inputs:
    %   ZIn   - Initial conditions (position and velocity)
    %   idx   - Index to start the forward and backward evolution
    %   dt    - Time step size
    %   tf    - Final simulation time
    %   eta   - Pre-generated stochastic terms
    %   Gamma - Damping coefficient
    %   kT    - Thermal energy (Boltzmann constant * Temperature)
    %   h4    - Quartic potential coefficient
    %   h2    - Quadratic potential coefficient
    % Outputs:
    %   outPath - Final trajectory path (position)
    %   outZ    - Final state of the system (position and velocity for each time step)

    % Time stepping parameters
    N = round(tf / dt); % Number of time steps

    %% Evolution Loop

    Sqdt = sqrt(dt); % Square root of time step for stochastic term
    Z = zeros(2, N); % Initialize trajectory
    Z(:, idx) = ZIn; % Set initial condition at given index

    % Forward evolution loop
    for n = idx + 1:N
        eta0 = eta(n-1);
        Z0 = Z(:, n-1);
        Z(:, n) = Z0 + DriftZ(Z0, Gamma, h4, h2) * dt + Stoch(Gamma, kT) * Sqdt * eta0;
    end

    % Backward evolution loop
    Z0 = ZIn;
    Z0(2) = -Z0(2); % Reverse velocity for backward evolution
    for n = idx - 1:-1:1
        eta0 = eta(n);
        Z0 = Z0 + DriftZ(Z0, Gamma, h4, h2) * dt + Stoch(Gamma, kT) * eta0 * Sqdt; % Euler method
        Z(:, n) = Z0;
    end

    % Reverse velocity for the backward segment
    Z(2, 1:idx-1) = -Z(2, 1:idx-1);

    %% Output
    outPath = Z(1, :); % Position trajectory
    outZ = Z; % Full state trajectory (position and velocity)
end

%% Drift Function for Quartic Well
function [out] = DriftZ(Z, Gamma, h4, h2)
    % Drift function for Langevin dynamics in a quartic potential well
    % Inputs:
    %   Z     - State vector [position; velocity]
    %   Gamma - Damping coefficient
    %   h4    - Quartic potential coefficient
    %   h2    - Quadratic potential coefficient
    % Output:
    %   out   - Drift vector
    out = [Z(2); 2 * h2 * Z(1) - 4 * h4 * Z(1)^3 - 2 * Gamma * Z(2)];
end

%% Stochastic Term Function
function [out] = Stoch(Gamma, Temp)
    % Stochastic function for Langevin dynamics
    % Inputs:
    %   Gamma - Damping coefficient
    %   Temp  - Temperature
    % Output:
    %   out   - Stochastic term vector
    out = 2 * sqrt(Gamma * Temp) * [0; 1];
end
