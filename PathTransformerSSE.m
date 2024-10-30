function [newPath, newPsi] = PathTransformerSSE(path, Psi, choice)
% Applies a specific transformation to the path and Hilbert-space trajectory based on the specified choice.
% Inputs:
%   path - Original path of the particle
%   Psi - Original Hilbert-space trajectory in Fock basis (bSize x trajectorySteps)
%   choice - Integer specifying the transformation type:
%       1 - No transformation
%       2 - Time reversal
%       3 - Parity (spatial reflection) reversal
%       4 - Parity-time (PT) reversal
% Outputs:
%   newPath - Transformed path
%   newPsi - Transformed Hilbert-space trajectory

    % Define parameters for parity transformation
    bSize = size(Psi, 1);
    trajectorySteps = size(Psi, 2);

    % Construct parity operator in Fock basis
    Parity = ((-ones(1, bSize)) .^ (0:bSize-1)).';
    Parity = repmat(Parity, 1, trajectorySteps); % Expand parity across trajectory steps

    % Apply transformation based on choice
    switch choice
        case 1
            % No transformation applied
            newPath = path;
            newPsi = Psi;
        case 2
            % Time reversal transformation
            [newPath, newPsi] = T_Reversal(path, Psi);
        case 3
            % Parity (P) reversal transformation
            [newPath, newPsi] = P_Reversal(path, Psi, Parity);
        case 4
            % Parity-time (PT) reversal transformation
            [newPath, newPsi] = PT_Reversal(path, Psi, Parity);
        otherwise
            error('Invalid choice: Please select 1, 2, 3, or 4.');
    end
end

%% Transformation Functions

function [outPath, outPsi] = T_Reversal(path, Psi)
% Applies time reversal to the Hilbert-space trajectory and path
% Inputs:
%   path - Original path of the particle
%   Psi - Original Hilbert-space trajectory in Fock basis
% Outputs:
%   outPath - Time-reversed path
%   outPsi - Time-reversed Hilbert-space trajectory

    outPsi = conj(flip(Psi, 2)); % Flip Psi along time axis and take complex conjugate
    outPath = flip(path); % Flip path along time axis
end

function [outPath, outPsi] = P_Reversal(path, Psi, Parity)
% Applies parity (spatial reflection) reversal to the Hilbert-space trajectory and path
% Inputs:
%   path - Original path of the particle
%   Psi - Original Hilbert-space trajectory in Fock basis
%   Parity - Parity operator matrix
% Outputs:
%   outPath - Parity-reversed path
%   outPsi - Parity-reversed Hilbert-space trajectory

    outPsi = Parity .* Psi; % Apply parity transformation to Psi
    outPath = -path; % Invert path (spatial reflection)
end

function [outPath, outPsi] = PT_Reversal(path, Psi, Parity)
% Applies parity-time (PT) reversal to the Hilbert-space trajectory and path
% Inputs:
%   path - Original path of the particle
%   Psi - Original Hilbert-space trajectory in Fock basis
%   Parity - Parity operator matrix
% Outputs:
%   outPath - PT-reversed path
%   outPsi - PT-reversed Hilbert-space trajectory

    outPsi = conj(flip(Parity .* Psi, 2)); % Apply parity, then flip along time, then take complex conjugate
    outPath = -flip(path); % Flip path along time and invert spatially
end
