function Pacc = PaccLangevinTIS(PsiNF, PsiOF, idx, idx_new, kT, Gamma, dt, h4, h2)
    % This function calculates the acceptance probability for Langevin dynamics in TIS.
    % Inputs:
    %   PsiNF, PsiOF - New and old path states
    %   idx, idx_new - Indices of the current and new shooting points
    %   kT           - Thermal energy (Boltzmann constant * Temperature)
    %   Gamma        - Damping coefficient
    %   dt           - Time step size
    %   h4, h2       - Potential coefficients for quartic and quadratic terms
    % Output:
    %   Pacc         - Acceptance probability

    % Determine the length of the new and old paths
    N_new = size(PsiNF, 2) - 2;
    N_old = size(PsiOF, 2) - 2;

    % Truncate the new and old paths to their respective indices
    PsiOF = PsiOF(:, 1:idx);
    PsiNF = PsiNF(:, 1:idx_new);

    % Define Hamiltonian function
    H = @(Z) (0.5 * Z(2)^2 + h4 * Z(1)^4 - h2 * Z(1)^2);

    % Initial acceptance probability based on energy difference
    P0 = exp(-(H(PsiNF(:, 1)) - H(PsiOF(:, 1))) / kT);

    % Time-reverse the new and old paths
    PsiNB = T_Reversal(PsiNF(:, 1:idx_new));
    PsiOB = T_Reversal(PsiOF(:, 1:idx));
    idx_max = max(idx, idx_new);

    % Initialize stochastic action variables
    SNF = zeros(1, idx_max);
    SNB = zeros(1, idx_max);
    SOF = zeros(1, idx_max);
    SOB = zeros(1, idx_max);
    sigma = [0; 2 * sqrt(Gamma * kT)];

    % Calculate stochastic action for new forward and backward paths
    for n = 1:idx_new - 1
        % New forward path
        Psi0 = PsiNF(:, n);
        mu = muFun(Psi0, Gamma, h4, h2);
        top = (sigma).' * ((PsiNF(:, n + 1) - Psi0) / dt - mu);
        SNF(n) = (top.' * top) / 2 / (sigma.' * sigma)^2;

        % New backward path
        Psi0 = PsiNB(:, n);
        mu = muFun(Psi0, Gamma, h4, h2);
        top = (sigma).' * ((PsiNB(:, n + 1) - Psi0) / dt - mu);
        SNB(n) = (top.' * top) / 2 / (sigma.' * sigma)^2;
    end

    % Calculate stochastic action for old forward and backward paths
    for n = 1:idx - 1
        % Old forward path
        Psi0 = PsiOF(:, n);
        mu = muFun(Psi0, Gamma, h4, h2);
        top = (sigma).' * ((PsiOF(:, n + 1) - Psi0) / dt - mu);
        SOF(n) = (top.' * top) / 2 / (sigma.' * sigma)^2;

        % Old backward path
        Psi0 = PsiOB(:, n);
        mu = muFun(Psi0, Gamma, h4, h2);
        top = (sigma).' * ((PsiOB(:, n + 1) - Psi0) / dt - mu);
        SOB(n) = (top.' * top) / 2 / (sigma.' * sigma)^2;
    end

    % Compute total stochastic action difference
    S = dt * sum(SNF - SNB + SOB - SOF);

    % Calculate the acceptance probability
    Pacc = real((N_old / N_new) * P0 * exp(-S));
end

% Time-Reversal Function
function [outPsi] = T_Reversal(Psi)
    % Reverses the trajectory and negates velocity
    Z = flip(Psi, 2);
    Z(2, :) = -Z(2, :);
    outPsi = Z;
end

% Drift Function for Quartic Well
function [outPsi] = muFun(Z, Gamma, h4, h2)
    % Drift function for Langevin dynamics in a quartic potential well
    % Inputs:
    %   Z     - State vector [position; velocity]
    %   Gamma - Damping coefficient
    %   h4    - Quartic potential coefficient
    %   h2    - Quadratic potential coefficient
    % Output:
    %   outPsi - Drift vector
    outPsi = [Z(2); 2 * h2 * Z(1) - 4 * h4 * Z(1)^3 - 2 * Gamma * Z(2)];
end
