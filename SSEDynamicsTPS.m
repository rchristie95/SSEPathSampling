function [outR, outPsi] = SSEDynamicsTPS(PsiIn, dt, tf, eta, Hhat, Lhat, Rhat, hbar, idx)
% Simulates the dynamics of a quantum state vector under the Stochastic Schrödinger Equation (SSE) with a given Hamiltonian, Lindblad operator, and reactive coordinate operator.
% Inputs:
%   PsiIn - Initial wavefunction (column vector)
%   dt - Time step size
%   tf - Final time of simulation
%   eta - Stochastic noise vector
%   Hhat - Hamiltonian operator
%   Lhat - Lindblad operator
%   Rhat - Reactive coordinate operator
%   hbar - Reduced Planck's constant
%   idx - Index separating backward and forward time evolution
% Outputs:
%   outR - Reactive coordinate at each time step
%   outPsi - Evolved wavefunction at each time step

    % Initialize variables
    bSize = size(PsiIn, 1); % Size of the state vector basis
    N = round(tf / dt); % Total number of time steps
    Sqdt = sqrt(dt); % Square root of time step for stochastic term
    Psi = zeros(bSize, N); % Matrix to store wavefunction at each time step
    ReactCoord = zeros(1, N); % Vector to store reactive coordinate

    % Initial conditions
    Psi(:, 1) = conj(PsiIn); % Initialize Psi with conjugated initial wavefunction
    ReactCoord(1) = real(PsiIn' * Rhat * PsiIn / (PsiIn' * PsiIn)); % Initial reactive coordinate
    Psi0 = Psi(:, 1); % Temporary variable for wavefunction update

    % Backward time evolution up to index 'idx'
    for n = 2:idx
        eta0 = eta(n-1); % Stochastic noise at previous step
        normPsi = norm(Psi0); % Normalize the wavefunction
        Psi0 = Psi0 + Drift(Psi0, Hhat, Lhat, bSize, hbar) * dt + Stoch(Psi0, Lhat, bSize, hbar) * eta0 * Sqdt; % Update wavefunction
        Psi0 = Psi0 / normPsi; % Normalize to unit norm
        Psi(:, n) = Psi0; % Store updated wavefunction
        ReactCoord(n) = real((Psi0' * Rhat * Psi0) / (Psi0' * Psi0)); % Update reactive coordinate
    end

    % Reverse the backward segment of the wavefunction and reactive coordinate
    Psi(:, 1:idx) = conj(flip(Psi(:, 1:idx), 2));
    ReactCoord(1:idx) = flip(ReactCoord(1:idx));

    % Forward time evolution from 'idx+1' to 'N'
    for n = idx+1:N
        eta0 = eta(n-1); % Stochastic noise at previous step
        Psi0 = Psi(:, n-1); % Previous wavefunction
        normPsi = norm(Psi0); % Normalize the wavefunction
        Psi0 = Psi0 + Drift(Psi0, Hhat, Lhat, bSize, hbar) * dt + Stoch(Psi0, Lhat, bSize, hbar) * eta0 * Sqdt; % Update wavefunction
        Psi0 = Psi0 / normPsi; % Normalize to unit norm
        Psi(:, n) = Psi0; % Store updated wavefunction
        ReactCoord(n) = real((Psi0' * Rhat * Psi0) / (Psi0' * Psi0)); % Update reactive coordinate
    end

    % Output final wavefunction and reactive coordinate
    outPsi = Psi; % Evolved wavefunction over time
    outR = real(ReactCoord); % Reactive coordinate over time
end

%% Supporting Functions

function [out] = Drift(Psi, H, L, bSize, hbar)
% Calculates the drift term for the wavefunction evolution
% Inputs:
%   Psi - Current wavefunction
%   H - Hamiltonian operator
%   L - Lindblad operator
%   bSize - Basis size of the wavefunction
%   hbar - Reduced Planck's constant
% Output:
%   out - Drift term for the SSE

    % Calculate expectation value of Lindblad operator
    EPsi = Psi' * L * Psi;
    % Compute the drift term
    out = hbar^(-1) * (-1i * H - 0.5 * (L' * L) + EPsi' * L - 0.5 * (EPsi * EPsi') * eye(bSize)) * Psi;
end

function [out] = Stoch(Psi, L, bSize, hbar)
% Calculates the stochastic term for the wavefunction evolution
% Inputs:
%   Psi - Current wavefunction
%   L - Lindblad operator
%   bSize - Basis size of the wavefunction
%   hbar - Reduced Planck's constant
% Output:
%   out - Stochastic term for the SSE

    % Compute the stochastic term
    out = (L - (Psi' * L * Psi) * eye(bSize)) * Psi * sqrt(1 / hbar);
end
