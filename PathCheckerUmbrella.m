function [choice, ns] = PathCheckerUmbrella(path, A, windowMat)
    % This function checks the validity of a given path in an umbrella sampling context and assigns a choice.
    % Inputs:
    %   path      - The given trajectory path
    %   A         - Boundaries defining state A
    %   windowMat - Matrix defining window boundaries for umbrella sampling
    % Outputs:
    %   choice    - Selected choice based on path validation (0 if invalid)
    %   ns        - Number of valid states (non-zero elements in Pdist)

    % Initialize probability distribution for path choices
    Pdist = [0, 0, 0, 0];
    trajectorySteps = length(path);

    % Check if the path is valid in different contexts
    if path(1) < A(2) && path(1) > A(1) && windowMat(1) < path(trajectorySteps) && path(trajectorySteps) < windowMat(2)
        Pdist(1) = 1; % Forward path within the window
    end

    if path(trajectorySteps) < A(2) && path(trajectorySteps) > A(1) && windowMat(1) < path(1) && path(1) < windowMat(2)
        Pdist(2) = 1; % Reversed path within the window
    end

    if path(1) > -A(2) && path(1) < -A(1) && -windowMat(1) > path(trajectorySteps) && path(trajectorySteps) > -windowMat(2)
        Pdist(3) = 1; % Mirrored forward path within the window
    end

    if path(trajectorySteps) > -A(2) && path(trajectorySteps) < -A(1) && -windowMat(1) > path(1) && path(1) > -windowMat(2)
        Pdist(4) = 1; % Mirrored reversed path within the window
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