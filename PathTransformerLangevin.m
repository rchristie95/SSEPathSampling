function [newPath, newPsi] = PathTransformerLangevin(path, Psi, choice)
% Applies a transformation to the input path and phase-space trajectory based on the specified choice
% Inputs:
%   path - Original path of the particle
%   Psi - Original phase-space trajectory of the particle
%   choice - Integer indicating the type of transformation:
%       1 - No transformation
%       2 - Time reversal
%       3 - Parity (spatial reflection) reversal
%       4 - Parity-time (PT) reversal
% Outputs:
%   newPath - Transformed path
%   newPsi - Transformed phase-space trajectory

    switch choice
        case 1
            % No transformation applied
            newPath = path;
            newPsi = Psi;
        case 2
            % Time reversal transformation
            [newPath, newPsi] = T_Reversal(Psi);
        case 3
            % Parity (P) reversal transformation
            [newPath, newPsi] = P_Reversal(Psi);
        case 4
            % Parity-time (PT) reversal transformation
            [newPath, newPsi] = PT_Reversal(Psi);
        otherwise
            error('Invalid choice: Please select 1, 2, 3, or 4.');
    end
end

function [outPath, outZ] = T_Reversal(Psi)
% Applies time reversal to the phase-space trajectory
% Inputs:
%   Psi - Original phase-space trajectory of the particle
% Outputs:
%   outPath - Transformed path (first row of the reversed phase-space trajectory)
%   outZ - Entire transformed phase-space trajectory after time reversal

    % Flip the phase-space trajectory along the time axis and invert momentum
    Z = flip(Psi, 2); % Flip across time
    Z(2, :) = -Z(2, :); % Invert momentum (second row)
    outZ = Z;
    outPath = outZ(1, :); % Extract position (first row) as the new path
end

function [outPath, outZ] = P_Reversal(Psi)
% Applies parity (spatial reflection) reversal to the phase-space trajectory
% Inputs:
%   Psi - Original phase-space trajectory of the particle
% Outputs:
%   outPath - Transformed path (first row of the reversed phase-space trajectory)
%   outZ - Entire transformed phase-space trajectory after parity reversal

    Z = -Psi; % Invert both position and momentum
    outZ = Z;
    outPath = outZ(1, :); % Extract position (first row) as the new path
end

function [outPath, outZ] = PT_Reversal(Psi)
% Applies parity-time (PT) reversal to the phase-space trajectory
% Inputs:
%   Psi - Original phase-space trajectory of the particle
% Outputs:
%   outPath - Transformed path (first row of the reversed phase-space trajectory)
%   outZ - Entire transformed phase-space trajectory after PT reversal

    % Flip phase-space trajectory along the time axis and invert position
    Z = flip(Psi, 2); % Flip across time
    Z(1, :) = -Z(1, :); % Invert position (first row)
    outZ = Z;
    outPath = outZ(1, :); % Extract position as the new path
end
