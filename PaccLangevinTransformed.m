function Pacc = PaccLangevinTransformed(SPsiN, PsiO, idx, kT, Gamma, dt, h4, h2, choice)
    % This function calculates the acceptance probability for Langevin dynamics
    % with a transformed path based on a given choice of symmetry.
    % Inputs:
    %   SPsiN, PsiO - New and old transformed path states
    %   idx         - Index of the shooting point
    %   kT          - Thermal energy (Boltzmann constant * Temperature)
    %   Gamma       - Damping coefficient
    %   dt          - Time step size
    %   h4, h2      - Potential coefficients for quartic and quadratic terms
    %   choice      - Choice of symmetry for path transformation (1: 'S', 2: 'T', 3: 'P', 4: 'PT')
    % Output:
    %   Pacc        - Acceptance probability

    % Determine the length of the path
    N = size(PsiO, 2);
    if choice == 2 || choice == 4
        idxPrime = N - (idx - 1);
    else
        idxPrime = idx;
    end

    % Select appropriate transformation function
    ChoiceString = {'S', 'T', 'P', 'PT'};
    functionName = strcat(ChoiceString{choice}, '_Reversal');
    PsiN = feval(functionName, SPsiN);
    PsiSO = feval(functionName, PsiO);

    % Hamiltonian function
    H = @(Z) (0.5 * Z(2)^2 + h4 * Z(1)^4 - h2 * Z(1)^2);

    % Initial acceptance probability based on energy difference
    P0 = exp((-1 / kT) * (H(SPsiN(:, 1)) - H(PsiO(:, 1))));

    % Time-reverse the transformed paths
    PsiNB = T_Reversal(PsiN(:, 1:idx));
    PsiSOB = T_Reversal(PsiSO(:, 1:idxPrime));

    % Initialize stochastic action variables
    S_SN = zeros(1, N);
    S_SO = zeros(1, N);
    S_N = zeros(1, N);
    S_O = zeros(1, N);
    sigma = 2 * sqrt(Gamma * kT) * [0; 1];

    % Calculate stochastic action for transformed and original paths
    for n = 1:N-1
        % Transformed new forward path
        Psi0 = SPsiN(:, n);
        mu = muFun(Psi0, h4, h2, Gamma);
        top = (sigma).' * ((SPsiN(:, n + 1) - Psi0) / dt - mu);
        S_SN(n) = (top.' * top) / 2 / (sigma.' * sigma)^2;

        % Original old forward path
        Psi0 = PsiO(:, n);
        mu = muFun(Psi0, h4, h2, Gamma);
        top = (sigma).' * ((PsiO(:, n + 1) - Psi0) / dt - mu);
        S_O(n) = (top.' * top) / 2 / (sigma.' * sigma)^2;
    end

    % Calculate stochastic action for time-reversed and transformed paths
    for n = 1:idx-1
        % Transformed new backward path
        Psi0 = PsiNB(:, n);
        mu = muFun(Psi0, h4, h2, Gamma);
        top = (sigma).' * ((PsiNB(:, n + 1) - Psi0) / dt - mu);
        S_N(n) = (top.' * top) / 2 / (sigma.' * sigma)^2;
    end

    for n = idx:N-1
        % Transformed new forward path
        Psi0 = PsiN(:, n);
        mu = muFun(Psi0, h4, h2, Gamma);
        top = (sigma).' * ((PsiN(:, n + 1) - Psi0) / dt - mu);
        S_N(n) = (top.' * top) / 2 / (sigma.' * sigma)^2;
    end

    for n = 1:idxPrime-1
        % Transformed old backward path
        Psi0 = PsiSOB(:, n);
        mu = muFun(Psi0, h4, h2, Gamma);
        top = (sigma).' * ((PsiSOB(:, n + 1) - Psi0) / dt - mu);
        S_SO(n) = (top.' * top) / 2 / (sigma.' * sigma)^2;
    end

    for n = idxPrime:N-1
        % Transformed old forward path
        Psi0 = PsiSO(:, n);
        mu = muFun(Psi0, h4, h2, Gamma);
        top = (sigma).' * ((PsiSO(:, n + 1) - Psi0) / dt - mu);
        S_SO(n) = (top.' * top) / 2 / (sigma.' * sigma)^2;
    end

    % Compute total stochastic action difference
    S = dt * sum(S_SN - S_O + S_SO - S_N);

    % Calculate the acceptance probability
    Pacc = real(P0 * exp(-S));
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

% Symmetry Transformation Functions
function [outZ] = S_Reversal(Psi)
    % Identity transformation (no change)
    outZ = Psi;
end

function [outZ] = T_Reversal(Psi)
    % Time reversal: flip trajectory and negate velocity
    Z = flip(Psi, 2);
    Z(2, :) = -Z(2, :);
    outZ = Z;
end

function [outZ] = P_Reversal(Psi)
    % Parity transformation: negate position and velocity
    Z = -Psi;
    outZ = Z;
end

function [outZ] = PT_Reversal(Psi)
    % Parity-time transformation: flip trajectory and negate position
    Z = flip(Psi, 2);
    Z(1, :) = -Z(1, :);
    outZ = Z;
end
