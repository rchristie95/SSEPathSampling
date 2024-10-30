function [choice, ns] = PathCheckerUnconstrained(path, A, B)
    % This function checks the validity of a given path in an unconstrained sampling context and assigns a choice.
    % Inputs:
    %   path - The given trajectory path
    %   A    - Boundaries defining state A
    %   B    - Boundaries defining state B
    % Outputs:
    %   choice - Selected choice based on path validation (0 if invalid)
    %   ns     - Number of valid states (non-zero elements in Pdist)

    % Initialize probability distribution for path choices
    Pdist = [0, 0, 0, 0];
    trajectorySteps = length(path);

    % Check if the path is valid in different contexts
    if path(1) < A(2) && path(1) > A(1) && any(path > B(1))
        Pdist(1) = 1; % Forward path leaving state A and reaching state B
    end

    if path(trajectorySteps) < A(2) && path(trajectorySteps) > A(1) && any(path > B(1))
        Pdist(2) = 1; % Reversed path leaving state A and reaching state B
    end

    if path(1) > B(1) && path(1) < B(2) && any(path < A(2))
        Pdist(3) = 1; % Parity-reversed path starting in state B and reaching state A
    end

    if path(trajectorySteps) > B(1) && path(trajectorySteps) < B(2) && any(path < A(2))
        Pdist(4) = 1; % Parity-time-reversed path starting in state B and reaching state A
    end

    % Count the number of valid states
    ns = length(nonzeros(Pdist));

    % Determine the choice
    if ns == 0
        choice = 0; % No valid choice
    else
        choice = randsample(1:4, 1, true, Pdist); % Randomly sample a valid choice based on Pdist
    end
end
