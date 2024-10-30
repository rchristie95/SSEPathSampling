function [out] = LindbladNull(Hhat, Lhat, hbar)
% Computes the null steady-state density matrix for the Lindblad master equation.
% Inputs:
%   Hhat - Hamiltonian matrix
%   Lhat - Lindblad operator matrix
%   hbar - Reduced Planck constant
% Output:
%   out - Steady-state density matrix (Hermitian, trace = 1)

% Initialize identity matrix for system size
bSize = size(Hhat, 1);
Id = eye(bSize);

% Construct the superoperator matrix M representing the Lindblad equation
M = (-1i / hbar) * (kron(Id, Hhat) - kron(Hhat.', Id)) + ...
    (1 / hbar) * kron(conj(Lhat), Lhat) - ...
    (1 / (2 * hbar)) * (kron(Id, Lhat' * Lhat) + kron((Lhat' * Lhat).', Id));

% Find the null space of M using the eigenvector corresponding to the smallest eigenvalue
[V, D] = eig(M);
D = diag(D);
[~, ind] = min(abs(D)); % Identify the smallest eigenvalue's index
Vt = V(:, ind); % Eigenvector associated with the smallest eigenvalue

% Reshape Vt to form the density matrix and normalize its trace to 1
Rho = reshape(Vt, [], bSize);
Rho = Rho / trace(Rho); % Ensures trace(Rho) = 1

% Ensure Rho is Hermitian by averaging with its conjugate transpose
out = 0.5 * (Rho + Rho');
end
