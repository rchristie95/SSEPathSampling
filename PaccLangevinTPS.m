function Pacc = PaccLangevinTPS(PsiNF, PsiOF, idx, kT, Gamma, dt, h4, h2)
    % This function calculates the acceptance probability for Langevin dynamics in TPS.
    % Inputs:
    %   PsiNF, PsiOF - New and old path states
    %   idx          - Index of the shooting point
    %   kT           - Thermal energy (Boltzmann constant * Temperature)
    %   Gamma        - Damping coefficient
    %   dt           - Time step size
    %   h4, h2       - Potential coefficients for quartic and quadratic terms
    % Output:
    %   Pacc         - Acceptance probability

    % Hamiltonian function
    H = @(Z) (0.5 * Z(2)^2 + h4 * Z(1)^4 - h2 * Z(1)^2);

    % Initial acceptance probability based on energy difference
    P0 = exp((-1 / kT) * (H(PsiNF(:, 1)) - H(PsiOF(:, 1))));

    % Time-reverse the new and old paths
    PsiNB = T_Reversal(PsiNF(:, 1:idx));
    PsiOB = T_Reversal(PsiOF(:, 1:idx));

    % Initialize stochastic action variables
    SNF = zeros(1, idx);
    SNB = zeros(1, idx);
    SOF = zeros(1, idx);
    SOB = zeros(1, idx);
    sigma = 2 * sqrt(Gamma * kT) * [0; 1];

    % Calculate stochastic action for new and old forward and backward paths
    for n = 1:idx-1
        % New forward path
        Psi0 = PsiNF(:, n);
        mu = muFun(Psi0, h4, h2, Gamma);
        top = (sigma).' * ((PsiNF(:, n + 1) - Psi0) / dt - mu);
        SNF(n) = (top.' * top) / 2 / (sigma.' * sigma)^2;

        % New backward path
        Psi0 = PsiNB(:, n);
        mu = muFun(Psi0, h4, h2, Gamma);
        top = (sigma).' * ((PsiNB(:, n + 1) - Psi0) / dt - mu);
        SNB(n) = (top.' * top) / 2 / (sigma.' * sigma)^2;

        % Old forward path
        Psi0 = PsiOF(:, n);
        mu = muFun(Psi0, h4, h2, Gamma);
        top = (sigma).' * ((PsiOF(:, n + 1) - Psi0) / dt - mu);
        SOF(n) = (top.' * top) / 2 / (sigma.' * sigma)^2;

        % Old backward path
        Psi0 = PsiOB(:, n);
        mu = muFun(Psi0, h4, h2, Gamma);
        top = (sigma).' * ((PsiOB(:, n + 1) - Psi0) / dt - mu);
        SOB(n) = (top.' * top) / 2 / (sigma.' * sigma)^2;
    end

    % Compute total stochastic action difference
    S = dt * sum(SNF - SNB + SOB - SOF);

    % Calculate the acceptance probability
    PT = exp(-S);
    Pacc = real(P0 * PT);
end

% Time-Reversal Function
function [outPsi] = T_Reversal(Psi)
    % Reverses the trajectory and negates velocity
    Z = flip(Psi, 2);
    Z(2, :) = -Z(2, :);
    outPsi = Z;
end

% Drift Function for Quartic Well
function [out] = muFun(Z, h4, h2, Gamma)
    % Drift function for Langevin dynamics in a quartic potential well
    % Inputs:
    %   Z     - State vector [position; velocity]
    %   h4    - Quartic potential coefficient
    %   h2    - Quadratic potential coefficient
    %   Gamma - Damping coefficient
    % Output:
    %   out   - Drift vector
    out = [Z(2); 2 * h2 * Z(1) - 4 * h4 * Z(1)^3 - 2 * Gamma * Z(2)];
end
