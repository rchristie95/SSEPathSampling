function [outPath, outPsi, idx_new, check, etaUsed] = SSEDynamicsTISLeft(PsiIn, dt, tf, Hhat, Lhat, Rhat, hbar, A)
    % This function simulates the dynamics of a quantum system based on the stochastic Schrödinger equation for TIS (left boundary).
    % Inputs:
    %   PsiIn  - Initial wavefunction (column vector)
    %   dt     - Time step size
    %   tf     - Final simulation time
    %   Hhat   - Hamiltonian operator
    %   Lhat   - Lindblad operator
    %   Rhat   - Observable operator
    %   hbar   - Reduced Planck's constant
    %   A      - Boundary for path evaluation
    % Outputs:
    %   outPath - Final trajectory of observable expectation values
    %   outPsi  - Final state of the system (wavefunction)
    %   idx_new - Index of the original input state in the trajectory
    %   check   - Indicator if the trajectory completed successfully
    %   etaUsed - The stochastic terms used during the simulation

    bSize = size(PsiIn, 1); % Size of the Hilbert space
    N = round(tf / dt); % Number of time steps
    Sqdt = sqrt(dt); % Square root of time step for stochastic term
    path = nan(1, N); % Path trajectory
    Psi = nan(bSize, N); % Forward trajectory
    Psib = nan(bSize, N); % Backward trajectory
    Sb = 0; % Backward success indicator
    check = 0; % Indicator if trajectory completes

    %% Backward Integration
    Psi0 = conj(PsiIn); % Set initial conditions for backward integration (complex conjugate)
    Psib(:, 1) = Psi0;
    path(1) = real((Psi0' * Rhat * Psi0) / (Psi0' * Psi0));
    etaUsed = nan(1, N); % Array to store stochastic values used

    for n = 2:N-1 % Backward evolution loop
        eta0 = randn;
        etaUsed(n-1) = eta0;
        normPsi = norm(Psi0);
        Psi0 = Psi0 + Drift(Psi0, Hhat, Lhat, bSize, hbar) * dt + Stoch(Psi0, Lhat, bSize, hbar) * eta0 * Sqdt;
        Psi0 = Psi0 / normPsi; % Euler method
        Psib(:, n) = Psi0;
        path(n) = real((Psi0' * Rhat * Psi0) / (Psi0' * Psi0));
        if path(n) >= A
            Sb = 1;
            break;
        end
    end

    % Flip the backward path to get the forward path
    path(1:n) = flip(path(1:n));
    Psi(:, 1:n) = conj(flip(Psib(:, 1:n), 2));
    idx_new = n;

    %% Forward Integration (if backward trajectory reached A)
    if Sb == 1
        for n = idx_new+1:N % Forward evolution loop
            eta0 = randn;
            etaUsed(n-1) = eta0;
            Psi0 = Psi(:, n-1);
            normPsi = norm(Psi0);
            Psi0 = Psi0 + Drift(Psi0, Hhat, Lhat, bSize, hbar) * dt + Stoch(Psi0, Lhat, bSize, hbar) * eta0 * Sqdt;
            Psi0 = Psi0 / normPsi; % Euler method
            Psi(:, n) = Psi0;
            path(n) = real((Psi0' * Rhat * Psi0) / (Psi0' * Psi0));
            if path(n) >= A
                check = 1;
                break;
            end
        end
    end

    % Remove unused elements from etaUsed
    etaUsed(isnan(etaUsed)) = [];

    %% Final Checks and Output Preparation
    if check == 1
        validData = ~isnan(path);
        Psi = Psi(:, validData); % Keep only the columns in Psi that correspond to non-NaN values in path
        path = path(validData); % Keep only the non-NaN values in path
    else
        Psi = nan(bSize, 1);
        path = nan;
        check = 0;
    end

    % Additional check for path length
    if length(path) < 3
        Psi = nan(bSize, 1);
        path = nan;
        check = 0;
    end

    %% Output
    if check == 1
        outPath = path;
        outPsi = Psi;
    else
        outPath = nan;
        outPsi = nan(bSize, 1);
    end
    idx_new = find(all(bsxfun(@eq, Psi, PsiIn), 1), 1, 'first');
end

%% Drift Function
function [out] = Drift(Psi, H, L, bSize, hbar)
    % Drift function for the stochastic Schrödinger equation
    % Inputs:
    %   Psi   - State vector (wavefunction)
    %   H     - Hamiltonian operator
    %   L     - Lindblad operator
    %   bSize - Size of the Hilbert space
    %   hbar  - Reduced Planck's constant
    % Output:
    %   out   - Drift vector
    EPsi = Psi' * L * Psi;
    out = hbar^(-1) * (-1i * H - 0.5 * (L' * L) + EPsi' * L - 0.5 * (EPsi * EPsi') * eye(bSize)) * Psi;
end

%% Stochastic Term Function
function [out] = Stoch(Psi, L, bSize, hbar)
    % Stochastic term function for the stochastic Schrödinger equation
    % Inputs:
    %   Psi   - State vector (wavefunction)
    %   L     - Lindblad operator
    %   bSize - Size of the Hilbert space
    %   hbar  - Reduced Planck's constant
    % Output:
    %   out   - Stochastic term vector
    out = (L - (Psi' * L * Psi) * eye(bSize)) * Psi * sqrt(1 / hbar);
end
