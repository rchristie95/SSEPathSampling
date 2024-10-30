function [outPath, outZ, idx_new, check, Tran, etaUsed] = LangevinDynamicsTISRight(ZIn, dt, tf, Gamma, Temp, h4, h2, A, B)
    % This function simulates the Langevin dynamics for a system transitioning between two states A and B.
    % Inputs:
    %   ZIn   - Initial conditions (position and velocity)
    %   dt    - Time step size
    %   tf    - Final simulation time
    %   Gamma - Damping coefficient
    %   Temp  - Temperature for stochastic term
    %   h4    - Quartic potential coefficient
    %   h2    - Quadratic potential coefficient
    %   A, B  - Boundaries of the states for evaluation
    % Outputs:
    %   outPath  - Final trajectory path (position)
    %   outZ     - Final state of the system (position and velocity)
    %   idx_new  - Index of the original input state in the trajectory
    %   check    - Indicator if the trajectory completed successfully
    %   Tran     - Indicator of transition between A and B
    %   etaUsed  - The stochastic terms used during the simulation

    %% Initial Data Setup

    % Time stepping parameters
    N = round(tf / dt); % Number of time steps
    etaUsed = nan(1, N); % Array to store stochastic values used
    eta = randn(1, N); % Pre-generated random stochastic terms

    % Initialize trajectory and auxiliary arrays
    Sqdt = sqrt(dt);
    Z = nan(2, N); % Forward trajectory
    Zb = nan(2, N); % Backward trajectory
    
    % Transition flags
    Sf = 0; % Forward success indicator
    Sb = 0; % Backward success indicator
    Tran = 0; % Transition indicator

    %% Backward Evolution Loop

    % Set initial conditions for backward integration
    Z0 = ZIn;
    Z0(2) = -Z0(2); % Reverse velocity for backward integration
    Zb(:, 1) = Z0;

    for n = 2:N-1 % Backward evolution loop
        eta0 = eta(n);
        etaUsed(n-1) = eta0;
        Z0 = Z0 + DriftZ(Z0, Gamma, h4, h2) * dt + Stoch(Gamma, Temp) * eta0 * Sqdt; % Euler method
        Zb(:, n) = Z0;

        % Check if trajectory reached boundary A or B
        if Zb(1, n) <= A
            Sb = 1;
            break;
        elseif Zb(1, n) >= B
            Sb = 0;
            break;
        end
    end

    % Flip and reverse velocity for forward evolution
    Z(:, 1:n) = flip(Zb(:, 1:n), 2);
    Z(2, :) = -Z(2, :);
    idx_new = n;

    %% Forward Evolution Loop (if backward trajectory reached A)

    if Sb == 1
        for n = idx_new + 1:N % Forward evolution loop
            eta0 = eta(n);
            etaUsed(n-1) = eta0;
            Z0 = Z(:, n-1);
            Z(:, n) = Z0 + DriftZ(Z0, Gamma, h4, h2) * dt + Stoch(Gamma, Temp) * Sqdt * eta0;

            % Check if trajectory reached boundary A or B
            if Z(1, n) <= A
                Sf = 1;
                break;
            elseif Z(1, n) >= B
                Sf = 1;
                Tran = 1;
                break;
            end
        end
    end

    %% Final Checks and Output Preparation

    if Sf == 1 || Sb == 1
        Z(:, any(isnan(Z))) = [];
        check = 1;
    else
        Z = nan(2, 1);
        check = 0;
    end

    % Remove unused elements from etaUsed
    etaUsed(isnan(etaUsed)) = [];

    % Final path and state
    path = Z(1, :);

    if length(path) < round(0.005 / dt)
        Z = nan(2, 1);
        check = 0;
        Tran = 0;
    end

    %% Output
    outPath = Z(1, :);
    outZ = Z;
    idx_new = find(all(bsxfun(@eq, Z, ZIn), 1), 1, 'first');
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