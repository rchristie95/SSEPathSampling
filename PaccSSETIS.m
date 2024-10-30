function [Pacc] = PaccSSETIS(PsiNF, PsiOF, idx, idx_new, ThermalState, Hhat, Lhat, hbar, dt)
    % This function calculates the acceptance probability for a quantum stochastic Schrödinger equation.
    % Inputs:
    %   PsiNF, PsiOF - New and old path states
    %   idx, idx_new - Indices of the shooting points for the old and new paths
    %   ThermalState - Thermal state matrix
    %   Hhat         - Hamiltonian operator
    %   Lhat         - Lindblad operator
    %   hbar         - Reduced Planck's constant
    %   dt           - Time step size
    % Output:
    %   Pacc         - Acceptance probability

    % Get the size of the Hilbert space and path lengths
    bSize = size(Hhat, 1);
    N_new = size(PsiNF, 2) - 2;
    N_old = size(PsiOF, 2) - 2;

    % Truncate paths to their respective shooting point indices
    PsiOF = PsiOF(:, 1:idx);
    PsiNF = PsiNF(:, 1:idx_new);
    PsiNB = conj(flip(PsiNF(:, 1:idx_new), 2));
    PsiOB = conj(flip(PsiOF(:, 1:idx), 2));

    % Initial acceptance probability based on thermal state
    P0 = real((PsiNF(:, 1)' * ThermalState * PsiNF(:, 1)) / (PsiOF(:, 1)' * ThermalState * PsiOF(:, 1)));

    % Initialize stochastic action variables
    idx_max = max(idx, idx_new);
    SNF = zeros(1, idx_max);
    SNB = zeros(1, idx_max);
    SOF = zeros(1, idx_max);
    SOB = zeros(1, idx_max);

    %% New Path Loop
    for n = 1:idx_new - 1
        % New forward path
        Psi0 = PsiNF(:, n) / norm(PsiNF(:, n));
        lev = PsiNF(:, n)' * Lhat * PsiNF(:, n) * eye(bSize);
        sigma = (1 / sqrt(hbar)) * (Lhat - lev) * Psi0;
        mu = hbar^(-1) * (-1i * Hhat - 0.5 * (Lhat' * Lhat) + lev' * Lhat - 0.5 * (lev' * lev)) * Psi0;
        top = (sigma)' * ((PsiNF(:, n + 1) - Psi0) / dt - mu);
        SNF(n) = (top' * top) / 2 / (sigma' * sigma)^2;

        % New backward path
        Psi0 = PsiNB(:, n) / norm(PsiNB(:, n));
        lev = PsiNB(:, n)' * Lhat * PsiNB(:, n) * eye(bSize);
        sigma = (1 / sqrt(hbar)) * (Lhat - lev) * Psi0;
        mu = hbar^(-1) * (-1i * Hhat - 0.5 * (Lhat' * Lhat) + lev' * Lhat - 0.5 * (lev' * lev)) * Psi0;
        top = (sigma)' * ((PsiNB(:, n + 1) - Psi0) / dt - mu);
        SNB(n) = (top' * top) / 2 / (sigma' * sigma)^2;
    end

    %% Old Path Loop
    for n = 1:idx - 1
        % Old forward path
        Psi0 = PsiOF(:, n) / norm(PsiOF(:, n));
        lev = PsiOF(:, n)' * Lhat * PsiOF(:, n) * eye(bSize);
        sigma = (1 / sqrt(hbar)) * (Lhat - lev) * Psi0;
        mu = hbar^(-1) * (-1i * Hhat - 0.5 * (Lhat' * Lhat) + lev' * Lhat - 0.5 * (lev' * lev)) * Psi0;
        top = (sigma)' * ((PsiOF(:, n + 1) - Psi0) / dt - mu);
        SOF(n) = (top' * top) / 2 / (sigma' * sigma)^2;

        % Old backward path
        Psi0 = PsiOB(:, n) / norm(PsiOB(:, n));
        lev = PsiOB(:, n)' * Lhat * PsiOB(:, n) * eye(bSize);
        sigma = (1 / sqrt(hbar)) * (Lhat - lev) * Psi0;
        mu = hbar^(-1) * (-1i * Hhat - 0.5 * (Lhat' * Lhat) + lev' * Lhat - 0.5 * (lev' * lev)) * Psi0;
        top = (sigma)' * ((PsiOB(:, n + 1) - Psi0) / dt - mu);
        SOB(n) = (top' * top) / 2 / (sigma' * sigma)^2;
    end

    % Compute total stochastic action difference
    S = real(dt * sum(SNF - SNB + SOB - SOF));

    % Calculate the acceptance probability
    Pacc = real((N_old / N_new) * P0 * exp(-S));
end
